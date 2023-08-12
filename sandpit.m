addpath("./Polar CODEC");
addpath("./MODEMs");
clear variables

EsN0Bob = 60;
EsN0Eve = 60;

iteration = 10;

codingRate = 'Dynamic';%Const

N = 1024;%Total length of the coded message
j = 50;%Number of scrambled bits to use
l = 4; %Specify a list size of 4addpath("./Polar CODEC");
addpath("./MODEMs");
clear variables

EsN0Bob = 60;
EsN0Eve = 60;

iteration = 10;

codingRate = 'Dynamic';%Const

N = 1024;%Total length of the coded message
j = 0;%Number of scrambled bits to use
l = 4; %Specify a list size of 4
%CRC bits
crc_len = 24;%24 for For 24A

b_errors = 0;
e_errors = 0;

count = 0;
%Jammer can be either: Disabled , Half Duplex, if its neither then it defaults to Full Duplex.
system = systemModel( N , j , l , crc_len , 1e-3 , codingRate , EsN0Bob , EsN0Eve, [ "Rayleigh" , "Non-Multipath" ] , 'Disabled' );

malicious = true;

SCL_Decoders = SCLdecoder( N , log2(N) + 1 , system.k , j , l , crc_len , system.bob.decoder.frozen_bits , system.bob.decoder.useful_bits , []);

%while ( count < iteration )

    [symbols, mbits , enc] = system.generateMessage();
    
    [rxBob,rxEve,gainsBob,gainsEve,jamGains] = system.wiretap_send(symbols);

    [bobMsg, eveMsg , llrs] = system.wiretap_demodulate( rxBob , rxEve, gainsBob , gainsEve , jamGains);

    SCL_Decoders(llrs);

    dec = zeros(1,N);

    loop = 1;
    done = 0;
    
    while(done == 0)
        [bit, done] = SC_Decoder(llrs);
        dec(loop) = bit;
        loop = loop + 1;
    end


    b_errors = b_errors + sum(bobMsg ~= mbits);

    e_errors = e_errors + sum(eveMsg ~= mbits); 


    count = count + 1;
%end

disp(['     Bob:    ' ' Eve:     ' 'Ratio']);
disp([ b_errors/(count * N)   e_errors/(count * N)   b_errors/e_errors]);


%CRC bits
crc_len = 24;%24 for For 24A

b_errors = 0;
e_errors = 0;

count = 0;
%Jammer can be either: Disabled , Half Duplex, if its neither then it defaults to Full Duplex.
system = systemModel( N , j , l , crc_len , 1e-3 , codingRate , EsN0Bob , EsN0Eve, [ "Rayleigh" , "Non-Multipath" ] , 'Disabled' );

malicious = true;

while ( count < iteration )

    [symbols, mbits] = system.generateMessage();
    
    [rxBob,rxEve,gainsBob,gainsEve,jamGains] = system.untrustedRelay_send(symbols,malicious);

    [bobMsg, eveMsg] = system.untrustedRelay_demodulate( rxBob , rxEve, gainsBob , gainsEve , jamGains ,  malicious);

    b_errors = b_errors + sum(bobMsg ~= mbits);

    e_errors = e_errors + sum(eveMsg ~= mbits); 


    count = count + 1;
end

disp(['     Bob:    ' ' Eve:     ' 'Ratio']);
disp([ b_errors/(count * N)   e_errors/(count * N)   b_errors/e_errors]);

