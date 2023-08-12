%{
Polar decoder class



%}
classdef polarDecoder < handle
%Public Methods
    methods (Access = public)
        function obj = polarDecoder( new_N , new_k , new_j, new_list_size, new_crc_len )
            obj.N = new_N;
            obj.k = new_k;
            obj.j = new_j;
            obj.n = obj.N - obj.k;
            obj.list_size = new_list_size;
            obj.crc_len = new_crc_len;
            
            obj.current_bit = 1;

            obj.reliability_sequence_full = load('./Polar CODEC/rel_seq_16384.mat').sequences;
            obj.reliability_sequence = obj.reliability_sequence_full( obj.reliability_sequence_full <= new_N ); 

            obj.frozen_bits = obj.reliability_sequence( 1 : obj.n );
            obj.useful_bits = obj.reliability_sequence( obj.n + 1 : obj.N );

        end
        
        function set_reliability_seq( obj , new_seq )
            
            obj.reliability_sequence_full = new_seq;
            obj.reliability_sequence = obj.reliability_sequence_full( obj.reliability_sequence_full <= new_N ); 

            obj.frozen_bits = obj.reliability_sequence( 1 : obj.n );
            obj.useful_bits = obj.reliability_sequence( obj.n + 1 : obj.N );
            
        end
        
        function update(obj, new_N , new_k , new_j , new_list_size, new_crc_len )

            obj.N = new_N;
            obj.k = new_k;
            obj.j = new_j;
            obj.n = obj.N - obj.k;
            obj.list_size = new_list_size;
            obj.crc_len = new_crc_len;

            obj.reliability_sequence = obj.reliability_sequence_full( obj.reliability_sequence_full <= new_N ); 

            obj.frozen_bits = obj.reliability_sequence( 1 : obj.n );
            obj.useful_bits = obj.reliability_sequence( obj.n + 1 : obj.N );
            
        end
        
        %This is the root of the tree
        function message = decode(obj,received , scrambled_bit_seed)


            scrambled_bit_pattern = obj.generate_scrambled_bits(scrambled_bit_seed);

            %Use the Successive Cancellation List decoder.
            %valid_messages = polar_code_slc_decode(received,obj.N,obj.list_size,obj.frozen_bits,obj.useful_bits);
            valid_messages = obj.polar_code_decode_scl( received , scrambled_bit_pattern , obj.N , obj.list_size);

            msg_bits = valid_messages(:,1:obj.k-obj.crc_len);
            
            computed_crc = nrCRCEncode(msg_bits','24C')';
            
            computed_crc = computed_crc(:,obj.k-obj.crc_len+1 :end);
            received_crc = valid_messages(:,obj.k-obj.crc_len+1 :end);
            
            comp = zeros(1,obj.list_size);
            valid = [];
            for idx = 1:1:obj.list_size
                comp(idx) = mean(computed_crc(idx,:)~=received_crc(idx,:));
                if( isequal(computed_crc(idx,:), received_crc(idx,:)) )
                    valid = [valid idx];%List the valid bits
                end   
            end
            if(isempty(valid))%If the valid bits are empty then we just guess one
                lowest = mink(comp,1);
                index = find(comp == lowest);
                message = msg_bits(index(1),:);
            else
                message = msg_bits(valid(1),:);
            end
        end

    end            
%{
*************** Internal functions needed for the decoder *****************
%}
    methods (Access = private)
        %Successive Canncellation list decoder
        %Credit in part to Andrew Thangaraj of the IIT Madras for
        %inspiration during development of the SCL decoder algorithm.
        %Cited material can be found at
        %https://github.com/tallamjr/iit-madras-5G-standard, <31 Dec 2020>
        %Date cited: November 2021
        function message = polar_code_decode_scl(obj, beliefs , scrambled_bit_pattern , N , list_size)
            %The maximum depth of the binary tree is log2(Number of bits in
            %codeword) offset by 1
            max_depth = log2(N) + 1;
            %make our LLRs and codewords matracies
            LLRs       = zeros(list_size, max_depth , N);%LLRs at each stage
            codewords  = zeros(list_size, max_depth , N);
            %make our path metrics
            pathMetric = Inf(list_size,1);
            pathMetric(1) = 0;
            %The current state of every node in the binary tree
            node_state = zeros(1,2*N-1);

            %Copy the incoming beliefs into the LLRs matrix ensuring each decoder
            %has the LLRs

            LLRs(:,1,:) = repmat(beliefs,list_size,1,1);


            node = 1;
            depth = 1; %Start at the root
            done = 0;

            while ( done == 0 )

                if ( depth == max_depth ) %We are at a leaf node
                    decision_metric = squeeze(LLRs(:,depth,node));
                    if(any(obj.frozen_bits == node))
                        codewords(:,depth , node ) = 0;%Make a decision, since this bit is frozen set it to 0
                        %If the bit is frozen, we update the pathmetric only if the
                        %decision metric was negative as this would mean
                        %the decoder made a counter intuitive decision
                        if (decision_metric < 0)
                            pathMetric = pathMetric + abs(decision_metric);
                        end

                    else
                        if(any(obj.scrambled_bits == node))
                            codewords(:,depth , node ) = scrambled_bit_pattern( find(obj.scrambled_bits == node) );

                        else
                            %Look at the decision metrics, if they are negative we make
                            %a decision that they are 1.
                            decision = decision_metric < 0;
                            %Now we made a decision we make the opposite decision, thus
                            %the tree size doubles
                            pathMetric_2 = [pathMetric; pathMetric+abs(decision_metric)];
                            %Then take the <list_size> number of decoders with the smallest path metric and remove the rest.
                            
                            [pathMetric, positions] = mink(pathMetric_2,list_size);
    
                            position_1 = positions > list_size;
                            positions(position_1) = positions(position_1) - list_size;
                            %Update the current bit decisions to be those of
                            %the remaining SC Decoders
                            decision = decision(positions);
                            %invert the bit decisions when they are the
                            %opposite of the decision matric
                            decision(position_1) = not( decision(position_1) );
                            %Remove all the decoders which have a path metric
                            %that is too high
                            LLRs = LLRs(positions,:,:);
    
                            codewords = codewords(positions,:,:);
                            %Assign the new decision to the codewords of the
                            %decoders
                            codewords(:,depth,node) = decision;
                        end


                    end
                    if(node == N)%If we are at the last node
                        done = 1;
                    end
                    depth = depth - 1;
                    node  = ceil(node/2);
                else
                    %If we are not at a leaf node, we are at an interior node
                    %If we are in an interior node, it can have one of 3
                    %possible states, either sending bits to the left,
                    %right, or returning bits to its invoking node
                    depth_range = 2^(max_depth - depth);%How many LLRs should the current depth be given
                    elements_to_use = 1 + (node - 1)*depth_range : (node)*depth_range;%What elements of the LLR matrix should we use
                    node_position = (2^(depth-1)-1) + node + 2;
                    switch node_state(node_position)
                       case 0%Left State
                           %The range when we half the current LLRs
                            use_LLRs = squeeze(LLRs(:,depth , elements_to_use));
                            %Update the node state to be right state
                            node_state(node_position) = 1;
                            %Go to the left child of this node first
                            depth = depth + 1;
                            node  = 2 * node - 1;
                            
                            LLRs(:,depth, elements_to_use( 1 : depth_range/2 )) = obj.polar_f( use_LLRs(:,1:depth_range/2), use_LLRs(:,1+depth_range/2:end));

                       case 1%Right State
                           %The range when we half the current LLRs
                            use_LLRs = squeeze(LLRs(:,depth, elements_to_use));

                            node_state(node_position) = 2;

                            depth = depth + 1;
                            node  = 2 * node;
                            bits = squeeze(codewords(:,depth,elements_to_use(1 : depth_range/2)));
                            LLRs(:,depth, elements_to_use(1+depth_range/2:end) ) = obj.polar_g(use_LLRs(:,1:depth_range/2),use_LLRs(:,1+depth_range/2:end), bits );

                       case 2%Returning State
                            use_bits = squeeze(codewords(:,depth + 1, elements_to_use ));

                            codewords(:,depth, elements_to_use) = obj.polar_XOR( use_bits(:,1:depth_range/2) , use_bits(:,1+depth_range/2:end) );

                            if(depth <= max_depth - 1)
                                depth = depth - 1;
                                node  = ceil(node/2);
                            end
                    end %switch-case

                end %if(depth == full_depth) ... else
            end %While done == 0

            message = squeeze(codewords(:,max_depth,obj.useful_bits));

        end
        %Polar Coding 'f' function, used to compute beliefs for the left
        %child of an interior node
        function result = polar_f( obj , belief_A , belief_B )
            result = ( sign(belief_A) .* sign(belief_B) .* min( abs(belief_A) , abs(belief_B) ) );
        end
        %Polar Coding 'g' function, used to compute beliefs for the right
        %child of an interior node
        function result = polar_g( obj, belief_A , belief_B , left_result )
            result = ((-1).^left_result).*belief_A + belief_B;
        end
        %Polar coding XOR function, used when a child returns its codeword
        %to its parent
        function result = polar_XOR( obj , bits_left , bits_right)
            result = [ xor(bits_left,bits_right) bits_right ];
        end 
        %Pass in the seed we need for generating the correct scrambled bit
        %sequence.
        function bits = generate_scrambled_bits(obj , seed)
            %Set the RNG to use the seed we pass in
            rng(seed);
            %Generate the scrambled bits
            bits = randi( [0 , 1] , [1 obj.j] );
        end

    end

    properties
        reliability_sequence_full;%Vectors of integers
        reliability_sequence;%Vectors of integers

        frozen_bits;
        useful_bits;
        scrambled_bits;

        current_bit;

        N;%Scalar, length of the current handle codeword
        k;%Scalar, length of the message bits, this includes the CRC
        j;
        n;

        list_size;%Scalar, length of the Successive Cancellation list
        crc_len;%Length of the CRC appended to the message bits
    end
end

