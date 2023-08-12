addpath("./Polar CODEC");
addpath("./MODEMs");
clear variables

iterations = 10000;

codingRate = 'Constant';%Const

N = 512;%Total length of the coded message
l = 4; %Specify a list size of 4

%What form of digital modulation to use

ModulationScheme = '256QAM';%Supports M-PSK and M-QAM

%M-QAM or M-PSK, the bits per symbol is log2(M), so BPSK M=2 QPSK M=4
%256QAM M = 256.

%CRC bits
crc_len = 24;%24 for For 24A

system = systemModel( N, l, crc_len, ModulationScheme, codingRate , 0 , 0 );

%Bob is decremented by this until 0
%Eve is incremented by this until it was the initial EbNo for bob
EbNoStep = 1;
EbNoStart = 60;

steps = 1;

symbolsUsed = system.get_symbol_order(system.bps,ModulationScheme);

NoIncs = EbNoStart / EbNoStep;

DeltaEbNo       = (2 * EbNoStart + 1 ) / 2;

SecrecyCapacity = zeros(1,length(DeltaEbNo));
EquvRate        = zeros(1,length(DeltaEbNo));


EbN0Bob = EbNoStart;
EbN0Eve = -EbNoStart;

while( EbN0Bob > -EbNoStart-1 )
    
    b_errors = 0;
    b_blkerr = 0;
    e_errors = 0;
    e_blkerr = 0;

    dcmcEve = 0;
    dcmcBob = 0;
    EqRate = 0;

    count  = 0;
    
    
    system.updateSNRs(EbN0Bob,EbN0Eve);

    while (count < iterations)

            [symbols, mbits] = system.generateMessage();

            [rxBob,rxEve,gainsBob,gainsEve] = system.send(symbols);

            %[bobMsg, eveMsg] = system.demodulate( rxBob , rxEve, gainsBob , gainsEve );



            dcmcEve = dcmcEve + system.get_dcmc(rxEve',gainsEve,system.eve.EbNo);
            dcmcBob = dcmcBob + system.get_dcmc(rxBob',gainsBob,system.bob.EbNo);


            EqRate = EqRate + system.get_equivocation_rate(rxEve,gainsEve,system.eve);


%{
            b_error = sum(bobMsg ~= mbits);
            if(b_error > 0)  
                b_errors = b_errors + b_error;
                b_blkerr = b_blkerr + 1;   
            end

            e_error = sum(eveMsg ~= mbits);
            if(e_error > 0)  
                e_errors = e_errors + e_error; 
            end 
%}

            count = count + 1;

    end
    
    DeltaEbNo(steps) = EbN0Bob - EbN0Eve;
    SecrecyCapacity(steps) = (dcmcBob - dcmcEve) / count;
    EquvRate(steps) = EqRate / count;
    
    EbN0Bob = EbN0Bob - 1;
    EbN0Eve = EbN0Eve + 1;
    
    steps = steps + 1;

end
%Mean of the SC and EQ Rate to find the average of throughout the
%simulations

figure(1)

plot(DeltaEbNo,SecrecyCapacity,'-s','markerSize',5);
xlabel('{\Delta}E_bN_o between Bob and Eve [dB]');
ylabel('Secrecy Capacity [bits/symbol]');

xlim([-EbNoStart  EbNoStart]);
set(gca, 'XDir','reverse')
title('Secrecy Capacity againsts SNR per bit difference between Bob and Eve');
grid on

figure(2)

plot(-60:1:60,EquvRate,'-s','markerSize',1);
xlabel('E_bN_o of Eve [dB]');
ylabel('Equivocation Rate [shannons]');

xlim([-EbNoStart  EbNoStart]);
title('Equivocation Rate againsts SNR per bit for Eve');
grid on
ylim([0 1.1]);