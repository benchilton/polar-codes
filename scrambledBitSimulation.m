addpath("./Polar CODEC");
addpath("./MODEMs");
clear variables

num_bit_err = 100;%Keep simulating until we reach this many bit errors

codingRate = 'Dynamic';%Const

N = 1024;%Total length of the coded message
l = 4; %Specify a list size of 4

%CRC bits
crc_len = 24;%24 for For 24A

j = 5;

len = 20;

Data = zeros(4,len);

parfor step = 1:len

    count = 0;

    snr = 2*step;

    b_errors = 0;
    e_errors = 0;

    EbN0Bob = snr;
    EbN0Eve = snr;
    
    system = systemModel( N , j , l , crc_len , 1e-3  , codingRate , EbN0Bob , EbN0Eve, ...
                     [ "Rayleigh" , "Non-Multipath" ] , 'Full Duplex' );
    
    while ( (b_errors < num_bit_err) || (e_errors < num_bit_err) )
    
        [symbols, mbits] = system.generateMessage();
                
        [rxBob,rxEve,gainsBob,gainsEve] = system.send(symbols);
    
        [bobMsg, eveMsg] = system.demodulate( rxBob , rxEve, gainsBob , gainsEve );
    
        b_errors = b_errors + sum(bobMsg ~= mbits);
    
        e_errors = e_errors + sum(eveMsg ~= mbits); 
    
        count = count + 1;
    end

    Data(:,step) = [ snr ; system.bps ; b_errors / (count * N); e_errors / (count * N)];

end

save( ['data_j_' num2str(j) '_.mat'] , 'Data' );
