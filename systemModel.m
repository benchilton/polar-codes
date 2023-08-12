%{
Class that contains the complete system model of a Guassian Wiretap channel
model which supports the conditions of a Rayleigh and AWGN channel model.

Includes:
    - A single authorised user, bob.
    - A single eavesdropper, eve


%}
classdef systemModel < handle
%Public Methods
    methods (Access = public)

        function obj = systemModel(N, j , l, crcType, modulationParam , codingRate, EsN0Bob, EsN0Eve, multipathChannel, jammerMode)
            obj.N = N;
            obj.listSize = l;
            obj.crc = crcType;

            obj.codingRateType = codingRate;

            %Load the SNR lookup curve fitting models from file
            obj.load_fitting_models();

            if( isstring( modulationParam ) )
                obj.ModulationScheme = modulationParam;
            else
                obj.ModulationScheme = upper ( obj.snrLookup( modulationParam , EsN0Bob ) ); 
            end

            [obj.k, obj.bps, obj.codingRate] = obj.get_coding_rate(obj.ModulationScheme, obj.codingRate, obj.N, obj.crc);

            obj.m = obj.k - obj.crc;
            obj.j = j;%How many scrambled bits we have in each codeword
            %Compute how many padding bits are to be needed
            obj.noPaddingBits = obj.bps * ceil( obj.N / obj.bps ) - obj.N;
            obj.paddingBits = zeros( 1 , obj.noPaddingBits );
            
            obj.encoder = polarEncoder(obj.N,obj.k);

            obj.set_modem_objects(obj.ModulationScheme);%Build the MODEM objects

            obj.bob    = obj.default_channelUser_settings( EsN0Bob , multipathChannel );
            obj.eve    = obj.default_channelUser_settings( EsN0Eve , multipathChannel );

            obj.setup_jammer( jammerMode , EsN0Bob , EsN0Eve , multipathChannel );

            obj.compute_frozen_bit_selection();

        end
        
        function updateSNRs(obj, bobEsNo, eveEsNo)
            
            obj.eve.EsNo = eveEsNo;
            obj.eve.EbNo = obj.eve.EbNo - ( 10*log10( obj.bps * obj.codingRate ) );
            obj.eve.channelNoise.EsNo = obj.eve.EsNo;

            obj.eve.noiseVar = 1 ./ 10^(obj.eve.EsNo/10);
            
            obj.eve.demodulator.Variance = obj.eve.noiseVar;

            obj.bob.EsNo = bobEsNo;
            obj.bob.EbNo = obj.bob.EbNo - 10*log10( obj.bps * obj.codingRate);
            obj.bob.noiseVar = 1 ./ 10^(obj.bob.EsNo/10);
            obj.bob.channelNoise.EsNo = obj.bob.EsNo;
            
            obj.bob.demodulator.Variance = obj.bob.noiseVar;
            
        end
        
%-------------------- Simulation methods --------------------

        function errors = berSim( obj , noiseCh , fadingCh , EsN0)

            [generated, message] = obj.generateMessage();

            [y, h] = fadingCh(generated);
            y = noiseCh(y);

            y_eq = lteEqualizeMMSE( y , h , 1./(10^EsN0) );

            s_dest = obj.bob.demodulator(y_eq);
            s_dest = s_dest( (obj.noPaddingBits + 1):end )';

            x_hat = obj.bob.decoder.decode( s_dest , obj.bob.scrambled_bits_seed );

            errors = sum(x_hat ~= message);

        end

        %Generates a random message modulates it and return it
        function [symbols , mbits, polarbits] = generateMessage(obj)
            mbits = randi( [0 , 1] , [1 obj.m] );
            %Append the 3GCC 5G standard 24 bit CRC to our message
            mbitsCrc = nrCRCEncode(mbits','24C')';

            [ scrambledBits , seedUsed ] = obj.generate_scrambled_bits( obj.j );

            obj.bob.scrambled_bits_seed = seedUsed;

            %Polar encode the message
            polarbits = obj.encoder.encode( mbitsCrc , scrambledBits );%Encode with a CRC
            %Digitally modulate the encoded bites using the desired modulation
            
            %Relic from the past, i decided to move away the NR modulate

            polarbits = [ obj.paddingBits polarbits ];%Add padding
            
            symbols = complex(obj.modulator(polarbits'));

        end
        
        %Sends generated symbols across the wiretap channel
        %The jammer will also attempt to jam the signal the eavesdropper receives
        %This function simulates transmission across a wiretap channel
        function [rxBob, rxEve, gainsBob, gainsEve, jamGains] = wiretap_send(obj, symbols)

            %We get in some symbols, we then pass the symbols through eve
            %and bobs channel.

            %For Bob, y = h * x + n
            [rxBob, gainsBob] = obj.bob.channelFading(symbols);
            rxBob = obj.bob.channelNoise(rxBob);

            %For Eve, y = h * x + n
            [rxEve, gainsEve] = obj.eve.channelFading(symbols);
            rxEve = obj.eve.channelNoise(rxEve);
            %Generate out jamming Signal
            obj.jammingSeed = rng; %Save the jammer Seed

            jammed = 0;

            if(obj.jammerLen ~= 0)
            
                jammingSignal = obj.modSymOrder (  randi([1 obj.symbCount] , size(rxBob) ) )';
                jamEve  = obj.jammer{2}.channelFading(jammingSignal);
                jammed = obj.jammer{2}.channelNoise(jamEve);

            end

            rxEve = jammed + rxEve;
            %If its a half duplex jammer then bob will receive the jamming
            %signal too.
            if(obj.jammerHDuplex == true)
                [jamBob , jamGains]  = obj.jammer{1}.channelFading(jammingSignal);
                jamBob = obj.jammer{1}.channelNoise(jamBob);
                rxBob = jamBob + rxBob;
            else
                jamGains = [];
            end

            gainsEve = sum(gainsEve,2);%Sum Columns
            gainsBob = sum(gainsBob,2);%Sum Columns

        end
        %Demodulates transmission across a wiretap channel.
        function [decBob, decEve, y_bob, y_eve] = wiretap_demodulate(obj, bobSymbols, eveSymbols , bobGains , eveGains , jamGains )
            
            %Before equalisation
            %If half duplex then bobs rec will be:
            % y = h*x_m + g*x_j + n
            %Assuming we know g and x_j we can remove the jamming
            if(obj.jammerHDuplex == true)
                rng(obj.jammingSeed); %Load the jammer Seed
                jammingSignal = obj.modSymOrder (  randi([1 obj.symbCount] , size(bobSymbols) ) )';
                bobSymbols = bobSymbols - jamGains.*jammingSignal;
            end

            %Use MMSE Equalisation algorithm to equalise the symbols

            eqBob = lteEqualizeMMSE(bobSymbols , bobGains , obj.bob.noiseVar );%Noise Variance must not be in dB
            
            eqEve = lteEqualizeMMSE(eveSymbols , eveGains , obj.eve.noiseVar );%Noise Variance must not be in dB

            y_bob = obj.bob.demodulator(eqBob);
            y_eve = obj.eve.demodulator(eqEve);
            
            y_bob = y_bob( (obj.noPaddingBits + 1):end )';
            y_eve = y_eve( (obj.noPaddingBits + 1):end )';

            decBob = obj.bob.decoder.decode( y_bob , obj.bob.scrambled_bits_seed );
            decEve = obj.eve.decoder.decode( y_eve , 'shuffle' );
        end
        
        %Sends generated symbols across the untrusted relay channel
        %The jammer will also attempt to jam the signal the untrusted relay receives
        function [rxDest, rxRelay, gainsDest, gainsRelay, jamGains] = untrustedRelay_send(obj, symbols, malicious)
            %We get in some symbols, we then pass the symbols through eve
            %and bobs channel.
            %In this model, Eve is an untrusted relay whereas Bob is the
            %intended destination.

            jammed = 0;
            %If the jammer is currently active then generate the jamming
            %signal
            jamGains = [];
            if(obj.jammerLen ~= 0)
            
                obj.jammingSeed = rng; %Save the jammer Seed
                jammingSignal = obj.modSymOrder (  randi([1 obj.symbCount] , size(symbols) ) )';
                [jamRelay , jamGains]  = obj.jammer{2}.channelFading(jammingSignal);
                jammed = obj.jammer{2}.channelNoise(jamRelay);

            end

            %Source sends the desired symbols to the relay using the
            %channel model
            [rxRelay, gainsRelay] = obj.eve.channelFading(symbols);
            rxRelay = obj.eve.channelNoise(rxRelay);
            %The cooperative jammer will send a jamming signal to the relay
            rxRelay = jammed + rxRelay;
            %If the relay is not malicious then it will equalise the
            %received symbols before sending them to the destination,
            %otherwise it will send them without and signal restoration
            %which makes it much harder for the destination to correct
            %decode the receives symbols
            toDest = rxRelay;
            if( malicious == false )
                toDest = lteEqualizeMMSE(toDest , gainsRelay , obj.eve.noiseVar );%Noise Variance must not be in dB
            end

            %For Destination, it will received the symbols from the relay
            [rxDest, gainsDest] = obj.bob.channelFading(toDest);
            rxDest = obj.bob.channelNoise(rxDest);

            %If its a half duplex jammer then bob will receive the jamming
            %signal too.

            %If multipath return the effective path gains for each symbol
            gainsRelay = sum(gainsRelay,2);%Sum Columns
            gainsDest = sum(gainsDest,2);%Sum Columns

        end
        %Demodulates transmission across an untrusted relay channel
        function [decDest, decRelay, y_dest, y_relay] = untrustedRelay_demodulate(obj, destSymbols, relaySymbols , destGains , relayGains , jamGains , malicious )
            
            %Before equalisation
            %If half duplex then dests rec will be:
            % y = h*x_m + g*x_j + n
            %Assuming we know g and x_j we can remove the jamming
            if(obj.jammerHDuplex == true)
                rng(obj.jammingSeed); %Load the jammer Seed
                jammingSignal = obj.modSymOrder (  randi([1 obj.symbCount] , size(destSymbols) ) )';

                eqDest = lteEqualizeMMSE(destSymbols , destGains , obj.bob.noiseVar );%Noise Variance must not be in dB

                eqDest = eqDest - jamGains.*jammingSignal;

                if(malicious == true)
                    eqDest = lteEqualizeMMSE(eqDest , relayGains , obj.eve.noiseVar );
                end

            else

                eqDest = lteEqualizeMMSE(destSymbols , destGains , obj.bob.noiseVar );%Noise Variance must not be in dB

                eqDest = lteEqualizeMMSE(eqDest , relayGains , obj.eve.noiseVar );

            end

            eqRelay = lteEqualizeMMSE(relaySymbols , relayGains , obj.eve.noiseVar );%Noise Variance must not be in dB

            y_dest = obj.bob.demodulator(eqDest);
            y_relay = obj.eve.demodulator(eqRelay);
            
            y_dest = y_dest( (obj.noPaddingBits + 1):end )';
            y_relay = y_relay( (obj.noPaddingBits + 1):end )';

            decDest = obj.bob.decoder.decode( y_dest , obj.bob.scrambled_bits_seed );
            decRelay = obj.eve.decoder.decode( y_relay , 'shuffle' );
        end
%-------------------- Generic Communications related methods --------------------
        %Builds modulator objects and sets the number of padding bits
        %needed
        function setModulation(obj,ModulationScheme)
            obj.ModulationScheme = upper(ModulationScheme);
            
            set_modem_objects(obj, ModulationScheme);

            obj.bob = obj.default_channelUser_settings(EbNoBob , multipathChannel);
            obj.eve = obj.default_channelUser_settings(EbNoEve , multipathChannel);
            
            %Compute how many padding bits are to be needed
            obj.noPaddingBits = obj.bps * ceil( obj.N / obj.bps ) - obj.N;
            obj.paddingBits = zeros( 1 , obj.noPaddingBits );
         
        end
        %Returns the coding rate and message bits per symbol based on
        %coding rate type, modulation scheme, etc..
        function [mbps, bps, codingRate] = get_coding_rate(obj, modScheme, codingMode, messageLen, crcLength)
            
            bps = obj.getBitsPerSymbol(modScheme);
            if(strcmp(codingMode,'Constant'))
                mbps = messageLen/2 - crcLength;
            else

                if(bps == 1)
                    mbps = messageLen/2 - crcLength;
                else
                    mbps = ceil( messageLen * ( (bps -1) / bps ) - crcLength );
                end

            end

            mbps = mbps + crcLength;
            
            codingRate = mbps / messageLen;%Coding Rate is k/N
            
        end

        %Plots the constellation diagram for the passed in symbols, and
        %optionally plot the actual constellation points, if not disregard
        function plot_constellation(obj,symbolsIn , ConstPoints, titleContents)
            symbolsIn = complex(symbolsIn);
            plot(symbolsIn,'b.');

            hold on
            
            if(~isempty(ConstPoints) )

                plot( complex(ConstPoints), 'r*' );
                
                if(~isempty(titleContents) )
                    title(titleContents);
                end
                
            else
                
                plot( complex(obj.modSymOrder) , 'r*' );
                
                if(~isempty(titleContents) )
                    title(titleContents);
                end
                
            end
            grid on
            hold off
            ylabel('Quadrature');
            xlabel('In-Phase');
            xlim([-1.5 1.5]);
            ylim([-1.5 1.5]);
            axis padded
        end
        
%-------------------- Security related methods --------------------

        %This function performs SNR Lookup, basically we look at BER curves
        %for the Rayleigh + AWGN channel for a desired BER and returns the
        %modulation scheme that best fits the channel EsN0.
        function [ modulationToUse , bps ] = snrLookup( obj , desiredBER, mainEsNo  )

            validModulations = [ "BPSK" "QPSK" "8QAM" "16QAM" "32QAM" "64QAM" "128QAM" "256QAM" ];
            len = length(validModulations);
            %Use the BER v SNR curves for the Rayleigh fading and AWGN
            %channel to determine which digital modulation scheme should be
            %used.
            %The decision from the BER cruves for which M-PSK/M-QAM is to be used is based on the main
            %channel
            %Allocate a vector for the EsN0 for each modulation scheme
            requiredEsN0 = Inf*ones(1, length(validModulations) );
            %Compute the required EsN0 values for each modulation scheme
            parfor modIter = 1:(len - 2)
                requiredEsN0( modIter ) = obj.findBER( obj.fittingModel{1,modIter} , desiredBER );
            end
            modOrder = zeros(1,len+1);
            modOrder(1:1:len) = requiredEsN0 < mainEsNo;

            bps = 1;
            %Find the highest order modulation which would give a BER low
            %enough
            
            while ( (modOrder(bps) == 1) && (modOrder(bps+1) == 1) && (bps < len) )
                bps = bps + 1;
            end
            modulationToUse = validModulations( bps );
        end

        function BPS = getBitsPerSymbol(obj, stringIn)
            if(strcmp(stringIn,'BPSK'))
                BPS = 1;
            else
                if(strcmp(stringIn,'QPSK'))
                    BPS = 2;
                else%If its not BPSK or QPSK, then we gotta do more work
                    stringIn = erase(stringIn,'-');
                    stringIn = erase(stringIn,'QAM');
                    stringIn = erase(stringIn,'PSK');
                    M = str2double(stringIn);
                    BPS = log2(M);    
                end
            end
        end
        %J is the amount of scrambled bits per message
        function [scrambledbits, seedUsed] = generate_scrambled_bits(obj, J)
            seedUsed = rng;

            scrambledbits = randi( [0 , 1] , [1 J] );
        end
        %{
        Builds the Jammer Objects to be either Half or full duplex.
        Half Duplex :
            - Jammer is independant from Bob
            - Assumption that both the Jammer and Bob are familiar with the
              jamming signal
        Full Duplex :
            - Bob is both a receiver and transmitter
        %}
        function setup_jammer( obj , jammerMode , EsN0Bob , EsN0Eve , multipathChannel )

            obj.jammerLen = obj.symbCount;
            if( matches(jammerMode,'Disabled') )
                obj.jammerLen = 0;
            end

            if( matches(jammerMode,'Half Duplex') )
                obj.jammerHDuplex = true;
            end

            %Jammer channel
            obj.jammer{1} = obj.default_channelUser_settings( EsN0Bob , multipathChannel );
            obj.jammer{2} = obj.default_channelUser_settings( EsN0Eve , multipathChannel );

        end

        function relSeq = monte_carlo_estimation(obj,N,channelUser,M)
            %Send a message of all 0's, decode and compute BER for each bit
            %channel
            n = log2(N);
            
            %Symbols never change as we just send 0s
            bits = zeros(N,1);
            syms = obj.modulator(bits);
            
            channelBERs = zeros(1,N);
            
            
            for interations = 1:M %Repeat for more accuracy 

                %[rec, gains] = channelUser.channelFading(syms);
                rec = channelUser.channelNoise( syms );
                
                %eq = lteEqualizeMMSE( rec , gains , channelUser.noiseVar );
                
                llrs = channelUser.demodulator( rec );

                channelBERs = channelBERs + obj.monteCarloSC(llrs, N , n);
                
            end

            [~,relSeq] = sort(channelBERs,'ascend');

        end
        
        function result = compute_bhattacharyya_parameter(obj,variance,rayleigh,scale)
            %obj.modSymOrder
            if(rayleigh == true)
                expectation = 1.*sqrt( 4/(4-pi)*log(4*scale) );
            else
                expectation = obj.modSymOrder;
            end
            
            A = 1/( sqrt(2*pi) * variance );
            
            func = @(y) A .* exp( -((y.^2 + expectation.^2)./(2*variance^2)) );

            result = integral(func,-Inf,Inf);

            
        end    
            
    end
        
%{
*************** Internal functions needed for the system *****************
%}
    %Private Methods
    methods (Access = private)     

        function compute_frozen_bit_selection(obj)

            sequences  = load('./Polar CODEC/rel_seq_16384.mat').sequences;

            bob_order = sequences(1,:);
            bob_order = bob_order( bob_order <= obj.N);

            message_slots   = bob_order( obj.N - obj.k + 1 : end );
            scrambled_slots = bob_order( obj.N - (obj.k + obj.j) + 1 : obj.N - obj.k );
            frozen_slots    = bob_order( 1 : obj.N - (obj.k + obj.j)  );

            obj.encoder.frozen_bits    = frozen_slots;
            obj.encoder.useful_bits    = message_slots;
            obj.encoder.scrambled_bits = scrambled_slots;

            obj.bob.decoder.frozen_bits    = frozen_slots;
            obj.bob.decoder.useful_bits    = message_slots;
            obj.bob.decoder.scrambled_bits = scrambled_slots;

            obj.eve.decoder.frozen_bits    = frozen_slots;
            obj.eve.decoder.useful_bits    = message_slots;
            obj.eve.decoder.scrambled_bits = scrambled_slots;

        end

        function channelUser = default_channelUser_settings(obj,EsNo,multipathChannel)
                    
                channelUser = channelReceiver();
                
                channelUser.EsNo = EsNo;
                channelUser.EbNo = EsNo - ( 10*log(obj.bps) + 10*log10(obj.codingRate * obj.bps) );
                channelUser.noiseVar = 1 ./ 10^(channelUser.EsNo/10);
    
                channelUser.decoder = polarDecoder( obj.N , obj.k , obj.j , obj.listSize , obj.crc );
                
                channelUser.channelNoise = comm.AWGNChannel('NoiseMethod',"Signal to noise ratio (Es/No)",'EsNo',channelUser.EsNo);
                                                        
                                                        
                
                if( matches( multipathChannel(1)  , "Rayleigh" , 'IgnoreCase' , true ) )

                    if( matches( multipathChannel(2) , 'Multipath' , 'IgnoreCase' , true ) )
                
                        channelUser.channelFading = comm.RayleighChannel( ...
                                                        'SampleRate',10e3, ...
                                                        'PathDelays',[0 1.5e-9], ...
                                                        'AveragePathGains',[2 1], ...
                                                        'NormalizePathGains',true, ...
                                                        'MaximumDopplerShift',30, ...
                                                        'DopplerSpectrum',{doppler('Gaussian',0.6),doppler('Flat')}, ...
                                                        'RandomStream','mt19937ar with seed', ...
                                                        'Seed',22, ...
                                                        'PathGainsOutputPort',true);
                    else
                    
                        channelUser.channelFading = comm.RayleighChannel( ...
                                                        'SampleRate',10e3, ...
                                                        'PathDelays',[0], ...
                                                        'AveragePathGains',[2], ...
                                                        'NormalizePathGains',true, ...
                                                        'MaximumDopplerShift',30, ...
                                                        'DopplerSpectrum',{doppler('Gaussian',0.6)}, ...
                                                        'RandomStream','mt19937ar with seed', ...
                                                        'Seed',22, ...
                                                        'PathGainsOutputPort',true);
                    
                    end

                end
                if( matches( multipathChannel(1)  , "Rician" , 'IgnoreCase' , true ) )

                    if( matches( multipathChannel(2) , 'Multipath' , 'IgnoreCase' , true ) )
                
                            channelUser.channelFading = comm.RicianChannel('SampleRate',1e6, ...
                                                            'PathDelays',[0.0 0.5 1.2]*1e-6, ...
                                                            'AveragePathGains',[0.1 0.5 0.2], ...
                                                            'KFactor',2.8, ...
                                                            'DirectPathDopplerShift',5.0, ...
                                                            'DirectPathInitialPhase',0.5, ...
                                                            'MaximumDopplerShift',50, ...
                                                            'DopplerSpectrum',doppler('Bell', 8), ...
                                                            'RandomStream','mt19937ar with seed', ...
                                                            'Seed',73, ...
                                                            'PathGainsOutputPort',true);
                    else
                    
                        channelUser.channelFading = comm.RicianChannel('SampleRate',1e6, ...
                                                            'PathDelays',0.0, ...
                                                            'AveragePathGains',0.1, ...
                                                            'KFactor',2.8, ...
                                                            'DirectPathDopplerShift',5.0, ...
                                                            'DirectPathInitialPhase',0.5, ...
                                                            'MaximumDopplerShift',50, ...
                                                            'DopplerSpectrum',doppler('Bell', 8), ...
                                                            'RandomStream','mt19937ar with seed', ...
                                                            'Seed',73, ...
                                                            'PathGainsOutputPort',true);
                    
                    end

                end

                channelUser.demodulator = myQAMDemodulator( obj.bps , obj.symbCount, channelUser.noiseVar );       
                channelUser.demodulator.Variance = channelUser.noiseVar;
                    
        end

        function set_modem_objects(obj, ModScheme)
            
            ModScheme = erase(ModScheme,"-");
            if( strcmp("BPSK",ModScheme) )
                obj.bps = 1;
                obj.symbCount = 2^obj.bps;
            end
            if( strcmp("QPSK",ModScheme) )
                obj.bps = 2;
                obj.symbCount = 2^obj.bps;
            end
            if( strcmp("8PSK",ModScheme) )
                obj.bps = 3;
                obj.symbCount = 2^obj.bps;
            end
            if( contains(ModScheme, "QAM" , 'IgnoreCase' ,true) )%The passed in modulation scheme is qam
                ModScheme = erase(ModScheme,"QAM");
                M = str2double(ModScheme);
                if( ~isnan(M) )
                    obj.symbCount = M;
                end
                obj.bps = log2(obj.symbCount);
                
                obj.symbCount = M;
            end

                obj.modSymOrder = qammod( (1:1:obj.symbCount) - 1 , obj.symbCount,'UnitAveragePower',true );
                obj.modulator   = myQAMModulator( obj.bps , obj.symbCount );
            
        end
        
        function secrecyCapacities = compute_modulation_SC(obj,modulations,iterations)
            
            numSupported = length(modulations);
            
            dcmcLoop = 0;
            count = 1;
            
            dcmcEve = 0;
            dcmcBob = 0;
            
            secrecyCapacities = zeros(1,numSupported);
            
            while(count < numSupported + 1)
                
                obj.setModulation(modulations(count));

                while(dcmcLoop < iterations)
               
                    symbols = obj.generateMessage();

                    [rxBob,rxEve,gainsBob,gainsEve] = obj.send(symbols);

                    dcmcEve = dcmcEve + obj.get_dcmc(rxEve',gainsEve,obj.eve.EbNo);
                    dcmcBob = dcmcBob + obj.get_dcmc(rxBob',gainsBob,obj.bob.EbNo);
                
                    dcmcLoop = dcmcLoop + 1;
                
                end
                dcmcLoop = 0;
                secrecyCapacities(count) = dcmcBob - dcmcEve;
                count = count + 1;
            end
            %We've summed up all the capacities hence now we just need to
            %average them
            secrecyCapacities = secrecyCapacities / iterations;
            
        end
        
        function channelErr = monteCarloSC(obj , llrsIn , N , n )
            
            channelErr = zeros(1,N);
            max_depth = n+1;
            
            LLRs      = zeros(n+1,N);
            states = zeros(n+1,N);
            LLRs(1,:) = llrsIn;
            
            node = 1;
            depth = 1; %Start at the root
            done = 0;

            while (done == 0)
                
                if(depth == max_depth )%At a leaf
                    %If the bit decoded was not 0, then an error occured
                    channelErr(node) = ( LLRs( max_depth , node ) < 0 );
                    
                    if(node == N)%Last Node
                        done = 1;
                    end
                    
                    depth = depth - 1;
                    node  = ceil(node/2);
                    
                else%We are at an interior node
                    depth_range = 2^(max_depth - depth);
                    elements_to_use = 1 + (node - 1)*depth_range : (node)*depth_range;
                    node_pos = (2^(depth-1)-1) + node + 2;
                    switch states(node_pos)
                       case 0%Left State
                           %The range when we half the current LLRs
                            use_LLRs = LLRs(depth , elements_to_use);
                            %Update the node state to be right state
                            states(node_pos) = 1;
                            %Go to the left child of this node first
                            depth = depth + 1;
                            node  = 2 * node - 1;

                            LLRs(depth, elements_to_use( 1 : depth_range/2 )) = ...
                            sign( use_LLRs(1:depth_range/2) ) .* sign( use_LLRs(1+depth_range/2:end) ) ...
                            .* min( abs( use_LLRs(1:depth_range/2) ) , abs(use_LLRs(1+depth_range/2:end)) );

                       case 1%Right State
                           %The range when we half the current LLRs
                            use_LLRs = LLRs(depth, elements_to_use);

                            states(node_pos) = 2;

                            depth = depth + 1;
                            node  = 2 * node;
                            LLRs(depth, elements_to_use(1+depth_range/2:end) ) = use_LLRs(1:depth_range/2) + use_LLRs( (1+depth_range/2 ):end);
 
                       case 2%Returning State
 
                            if(depth <= max_depth - 1)
                                depth = depth - 1;
                                node  = ceil(node/2);
                            end
                    end %switch-case
                end
            end   
        end
            
        function load_fitting_models(obj)
            %There are 8 supported modulation schemes and hence 8 fitting models needed
            obj.fittingModel = cell(1,8);
            %The models are ordered as follows:
            %1:BPSK
            %2:QPSK
            %3:8PSK
            %4:16QAM
            %5:32QAM
            %6:64QAM
            %7:128QAM
            %8:256QAM
            %Placeholders
            obj.fittingModel{1,1} = load('./Fitting Models/model_bpsk.mat').fittedmodel;
            obj.fittingModel{1,2} = load('./Fitting Models/model_qpsk.mat').fittedmodel;
            obj.fittingModel{1,3} = load('./Fitting Models/model_8qam.mat').fittedmodel;
            obj.fittingModel{1,4} = load('./Fitting Models/model_16qam.mat').fittedmodel;
            obj.fittingModel{1,5} = load('./Fitting Models/model_32qam.mat').fittedmodel;
            obj.fittingModel{1,6} = load('./Fitting Models/model_64qam.mat').fittedmodel;
            obj.fittingModel{1,7} = load('./Fitting Models/model_128qam.mat').fittedmodel;
            obj.fittingModel{1,8} = load('./Fitting Models/model_256qam.mat').fittedmodel;

        end

        function EsN0 = findBER(obj , fitting , desiredBER )
            
            f = @(x) ( desiredBER - fitting(x) );
            EsN0 = fzero( f, [0 60]);

        end

    end

    properties
        
%---------- Doubles ----------

    N;%Total length of the polar code
    m;%Total length of the message bits excluding the crc
    k;%Total length of the message bits, this includes a crc
    j;%Number of scrambled bits
    
    listSize;%how big should the list be for the SCL decoder
    
    crc;%Length of the CRC
    
    
    bps;%bits per symbol
    symbCount;%2^bps is the number of symbols in the constellation
    
    modSymOrder;%Order of the symbols
    
    %How Many Padding Bit we need 
    noPaddingBits

    jammerHDuplex;
    jammerLen;
    jammingSeed;
%---------- Strings ----------

    ModulationScheme;%What modulation scheme this system uses
    codingRate;
    codingRateType;%Constant or Dynamic
    
%---------- Vectors ----------

    paddingBits;%The padding bits vector
    
%---------- Functions ----------
%---------- Classes ----------

    %----- Modulator -----
        %The decoders and demodulators will change however
        modulator;%This is either an MPSK or an MQAM Modulator

    %----- Send Encoder -----

        encoder;%Alice's Encoder on the wiretap channel
    
    %----- Users -----

        bob; %contains a decoder and channel models
        eve; %contains a decoder and channel models

        jammer;

        %Eve can easily be made into an array and hence multiple
        %eavesdroppers
    %----- Curve-fitting Models -----

        fittingModel;
        
    end
end

