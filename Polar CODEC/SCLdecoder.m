classdef SCLdecoder < matlab.System
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

        sc_decoders;

    end

    properties(DiscreteState)

    end

    % Pre-computed constants
    properties(Access = private)

    end
    
    methods(Access = public)
        
        function obj = SCLdecoder( N , n , k , j , l , crc , fb , ub , sb)

            obj.N = N;
            obj.n = n;%n = log2(N) - 1 to work
            obj.k = k;
            obj.j = j;
            obj.l = l;
            obj.crc = crc;
            obj.fb = fb;
            obj.ub = ub;
            obj.sb = sb;

            obj.sc_decoders = cell( 1 , 2*obj.l );

            for iter = 1:1:obj.l
                obj.sc_decoders{iter} = SCdecoder( obj.N , obj.n , obj.k , obj.j , obj.l , obj.crc , obj.fb , obj.ub , obj.sb );
            end

        end

    end

    methods(Access = protected)
        %Performs SC decoder operation until a leaf node makes a decision
        %The state is stored
      
        function decoded = stepImpl( obj , LLRs )

            decoding = true;

            decisions = zeros( 1,2 * obj.l );

            llrs = zeros( 1,2 * obj.l );

            path_metric = zeros( 1,2 * obj.l );

            while(decoding == true)

                %Each SC decoder performs SC decoding until a leaf node
                %makes a decision.
                for iter = 1:1:obj.l
                    decoder = obj.sc_decoders{iter};
                    [ decisions(iter) , llrs(iter) , ~ ] = decoder(LLRs)
                    obj.sc_decoders{obj.l + iter} = clone(decoder);
                end

                obj.sc_decoders = {obj.sc_decoders clone(obj.sc_decoders)}

                %Make the inverse bit decision of each decoder
                decisions(1+ end/2 : end) = not( decisions(1 : end/2 ) );

            end

        end
        
    end

    methods(Access = private)

    end

end
