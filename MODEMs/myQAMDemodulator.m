classdef myQAMDemodulator < matlab.System
    % Public, tunable properties
    properties
        BitsPerSymbol
        NoSymbols
        Variance
    end

    properties(DiscreteState)

    end

    % Pre-computed constants
    properties(Access = private)

    end
    
    methods(Access = public)
        
        function obj = myQAMDemodulator(bps,symbCount,NVar)
            obj.BitsPerSymbol = bps;
            obj.NoSymbols     = symbCount;
            obj.Variance      = NVar;%Not in dB
        end
        
    end

    methods(Access = protected)

        function LLRsOut = stepImpl(obj,symbsIn)
           LLRsOut = qamdemod(symbsIn,obj.NoSymbols,'OutputType','approxllr','UnitAveragePower',true,'NoiseVariance',obj.Variance);
        end
        
    end
end
