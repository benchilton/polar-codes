classdef myQAMModulator < matlab.System
    % Public, tunable properties
    properties
        BitsPerSymbol
        NoSymbols
    end

    properties(DiscreteState)

    end

    % Pre-computed constants
    properties(Access = private)

    end
    
    methods(Access = public)
        
        function obj = myQAMModulator(bps,symbCount)
            obj.BitsPerSymbol = bps;
            obj.NoSymbols     = symbCount;
        end
        
    end

    methods(Access = protected)

        function symbolsOut = stepImpl(obj,bitsIn)
           symbolsOut = qammod(bitsIn, obj.NoSymbols , 'InputType' , 'bit', 'UnitAveragePower',true );
        end
        
    end
end
