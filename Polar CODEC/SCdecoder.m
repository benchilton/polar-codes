classdef SCdecoder < matlab.System
    % Public, tunable properties
    properties
        N;
        n;
        k;
        j;
        l;
        crc;
        fb;
        ub;
        sb;

        node;
        depth;

        llrs;
        nodeStates;
        codes;
        clean_slate;

    end

    properties(DiscreteState)

    end

    % Pre-computed constants
    properties(Access = private)

    end
    
    methods(Access = public)
        
        function obj = SCdecoder( N , n , k , j , l , crc , fb , ub , sb)

            obj.N = N;
            obj.n = n;%n = log2(N) - 1 to work
            obj.k = k;
            obj.j = j;
            obj.l = l;
            obj.crc = crc;
            obj.fb = fb;
            obj.ub = ub;
            obj.sb = sb;

            obj.llrs       = zeros( n, N );
            obj.nodeStates = zeros(1,2*N-1);
            obj.codes = zeros( n, N);
            obj.node  = 0;
            obj.depth = 0;

            obj.clean_slate = 1;

        end

    end

    methods(Access = protected)
        %Performs SC decoder operation until a leaf node makes a decision
        %The state is stored
      
        function [next_bit , decision_metric , finished ] = stepImpl( obj , LLRs )
            
            next_decision = 0;

            if(obj.clean_slate == 1)
               %reset decoder state
               obj.llrs       = zeros( obj.n, obj.N );
               obj.nodeStates = zeros( obj.n, obj.N );
               obj.codes      = zeros( obj.n, obj.N);
               obj.node  = 1;
               obj.depth = 1;
               obj.llrs(1,:) = LLRs;
               obj.clean_slate = 0;
            end

            while( next_decision == 0 )
                if( obj.depth == obj.n )%leaf node

                    decision_LLR = obj.llrs(obj.depth,obj.node);

                    if(any(obj.fb == obj.node))%If the current bit position we are in is frozen
                        next_bit = 0;

                        decision_metric = abs(decision_LLR) * (decision_LLR < 0);

                    else
                        if( obj.llrs(obj.depth,obj.node) >= 0)
                           next_bit = 0;
                        else
                           next_bit = 1;
                        end
                        decision_metric = abs(decision_LLR) * (decision_LLR < 0);
                    end
                    finished = 0;
                    if( obj.node == obj.N )
                        finished = 1;
                        obj.clean_slate = 1;
                    end
                    %The decoder made a bit decision so pause
                    next_decision = 1;

                    obj.depth = obj.depth - 1;
                    obj.node  = ceil(obj.node/2);

                else
                    %If we are not at a leaf node, we are at an interior node
                    %If we are in an interior node, it can have one of 3
                    %possible states, either sending bits to the left,
                    %right, or returning bits to its invoking node
                    depth_range = 2^(obj.n - obj.depth);%How many LLRs should the current depth be given
                    elements_to_use = 1 + (obj.node - 1)*depth_range : (obj.node)*depth_range;%What elements of the LLR matrix should we use
                    node_position = (2^(obj.depth-1)-1) + obj.node + 2;
                    switch obj.nodeStates(obj.depth,obj.node)
                        case 0%Left State
                           %The range when we half the current LLRs
                            use_LLRs = obj.llrs(obj.depth , elements_to_use);
                            %Update the node state to be right state
                            obj.nodeStates(obj.depth,obj.node) = 1;
                            %Go to the left child of this node first
                            obj.depth = obj.depth + 1;
                            obj.node  = 2 * obj.node - 1;
                            
                            obj.llrs(obj.depth, elements_to_use( 1 : depth_range/2 )) = obj.polar_f( use_LLRs(:,1:depth_range/2), use_LLRs(:,1+depth_range/2:end));

                       case 1%Right State
                           %The range when we half the current LLRs
                            use_LLRs = obj.llrs(obj.depth, elements_to_use);

                            obj.nodeStates(obj.depth,obj.node) = 2;

                            obj.depth = obj.depth + 1;
                            obj.node  = 2 * obj.node;
                            bits = obj.codes(obj.depth,elements_to_use(1 : depth_range/2));
                            obj.llrs(obj.depth, elements_to_use(1+depth_range/2:end) ) = obj.polar_g(use_LLRs(:,1:depth_range/2),use_LLRs(:,1+depth_range/2:end), bits );

                       case 2%Returning State
                            use_bits = obj.codes(obj.depth + 1, elements_to_use );

                            obj.codes(obj.depth, elements_to_use) = obj.polar_XOR( use_bits(1:depth_range/2) , use_bits(1+depth_range/2:end) );

                            if(obj.depth <= obj.n - 1)
                                obj.depth = obj.depth - 1;
                                obj.node  = ceil(obj.node/2);
                            end

                    end %switch-case

                end

            end

        end
        
    end

    methods(Access = private)

        %Polar Coding 'f' function, used to compute beliefs for the left
        %child of an interior node
        function result = polar_f( ~ , belief_A , belief_B )
            result = ( sign(belief_A) .* sign(belief_B) .* min( abs(belief_A) , abs(belief_B) ) );
        end
        %Polar Coding 'g' function, used to compute beliefs for the right
        %child of an interior node
        function result = polar_g( ~, belief_A , belief_B , left_result )
            result = ((-1).^left_result).*belief_A + belief_B;
        end
        %Polar coding XOR function, used when a child returns its codeword
        %to its parent
        function result = polar_XOR( ~ , bits_left , bits_right)
            result = [ xor(bits_left,bits_right) bits_right ];
        end 

    end

end
