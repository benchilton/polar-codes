%{
Channel User Class

Contains different channel models

%}
classdef channelReceiver < handle
%Public Methods
    methods (Access = public)
        
        function obj = channelReceiver()

            obj.blockErrors = 0;
            obj.errors = 0;
            
        end

    end            
%{
*************** Internal functions needed for the decoder *****************
%}
    methods (Access = private)
 
    end
    
    properties

        scrambled_bits_seed;%Seed for the scrambled  bits
        
        demodulator;%This is either an MPSK or an MQAM Demodulator
        
        %Channel Models
        channelNoise;
        channelFading;

        %Polar Decoder:
        decoder;
        
        errors
        blockErrors
        
        EbNo;%Bit Power to Noise Ratio, in [dB]
        EsNo;%Symbol Power to Noise Ratio, in [dB]
        noiseVar;%Measured in regular units
        
    end
end

