%{
Polar encoder class



%}
classdef polarEncoder < handle
%Public Methods
    methods (Access = public)
        function obj = polarEncoder( new_N , new_k )

            obj.N = new_N;
            obj.k = new_k;
            obj.n = obj.N - obj.k;
            
            %obj.reliability_sequence_full = load('./Polar CODEC/rel_seq.mat').relSeq;
            obj.reliability_sequence_full = load('./Polar CODEC/rel_seq_16384.mat').sequences(1,:);
            obj.reliability_sequence = obj.reliability_sequence_full( obj.reliability_sequence_full <= new_N ); 

            obj.frozen_bits = obj.reliability_sequence( 1 : obj.n );
            obj.useful_bits = obj.reliability_sequence( obj.n + 1 : obj.N );
            
            
            
            obj.G_1 = [1 0 ; 1 1];
            obj.G_N = obj.G_1;
            for i = 1:1:(log2( obj.N ) - 1)
                obj.G_N = kron(obj.G_1 , obj.G_N);
            end
            
            
            
        end
        
        function obj = update( new_N , new_k )
            
            obj.N = new_N;
            obj.k = new_k;
            obj.n = obj.N - obj.k;
            
            obj.reliability_sequence = obj.reliability_sequence_full( obj.reliability_sequence_full <= new_N );
            obj.frozen_bits = obj.reliability_sequence( 1 : obj.n );
            obj.useful_bits = obj.reliability_sequence( obj.n + 1 : obj.N );
            
            
            %Make the generator matrix
            obj.G_N = obj.G_1;

            for i = 1:1:(log2( obj.N ) - 1)
                obj.G_N = kron(obj.G_1 , obj.G_N);
            end
            
        end
        
        %This is the root of the tree
        function codeword = encode(obj,message_bits , scrambled_sequence)
            
            codeword = zeros(1,obj.N);

            codeword(obj.scrambled_bits) = scrambled_sequence;

            codeword(obj.useful_bits) = message_bits;
            codeword = codeword * obj.G_N;
            codeword = mod(codeword,2);
            
        end

    end            
%{
*************** Internal properties needed for the encoder *****************
%}
    properties
        
        G_1; % [1 0 ; 1 1]
        G_N; %Generator matrix for the polar transform
        
        reliability_sequence_full;%Vectors of integers
        reliability_sequence;%Vectors of integers

        frozen_bits;
        useful_bits;
        scrambled_bits;

        N;%Scalar, length of the current handle codeword
        n;%Number of message bits
        k;%Scalar, length of the message bits, this includes the CRC

    end
end

