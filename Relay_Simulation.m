function Relay_Simulation( iterations , j_in , filename)

    name_to_use = ['data_' filename '_' num2str(j_in) '.csv'];
    
    if exist(name_to_use, 'file')==2
      delete(name_to_use);
    end

    codingRate = 'Dynamic';%Const
    
    N = 1024;%Total length of the coded message
    l = 4; %Specify a list size of 4
    
    %CRC bits
    crc_len = 24;%24 for For 24A
    
    %What form of digital modulation to use
    
    %M-QAM or M-PSK, the bits per symbol is log2(M), so BPSK M=2 QPSK M=4
    %256QAM M = 256.
    
    EsN0_start = 0;%EbNo Value to begin with
    EsN0_stop  = 60;%EbNo Value to end with
    EsN0_step  = 1.5;%Step for EbNo

    parfor index = EsN0_start:((EsN0_stop)/EsN0_step)

        EsN0 = index * EsN0_step;

        system = systemModel( N , j_in , l , crc_len , 1e-3 , codingRate , EsN0 , EsN0 , [ "Rayleigh" , "Non-Multipath" ] , 'Half Duplex' );
        
        b_errors = 0;
        e_errors = 0;
        
        count = 0;

        malicous = true;
        
    while ( (b_errors < iterations) || (e_errors < iterations) )
            
        [symbols, mbits] = system.generateMessage();
        
        [rxBob,rxEve,gainsBob,gainsEve,jamGains] = system.untrustedRelay_send(symbols,malicous);
    
        [bobMsg, eveMsg] = system.untrustedRelay_demodulate( rxBob , rxEve, gainsBob , gainsEve , jamGains ,  malicous);
    
        b_errors = b_errors + sum(bobMsg ~= mbits);

        e_errors = e_errors + sum(eveMsg ~= mbits); 

        count = count + 1;
    end
    
    to_write = [ ...
                 b_errors/(count * N) , e_errors/(count * N) , EsN0, ...
                 system.bps
                ];

    disp(to_write);

    writematrix(to_write, name_to_use ,'WriteMode','append');
    
    end
end
