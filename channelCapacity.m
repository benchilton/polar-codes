clear variables;
SNR_dB = -10:2:50;

% Affects the accuracy and duration of the simulation
symbol_count = 1e6;

power = 1; %16 is 16QAM, 256 is 256QAM, 2 is BPSK etc..
lim = 8;

% SNR range

snr_dB = SNR_dB;
snr = 10.^(snr_dB/10);
C = log2(1+snr);

% Modulation scheme

rayleigh_channel = comm.RayleighChannel( ...
    'SampleRate',10e3, ...
    'PathDelays',[0], ...
    'AveragePathGains',[2], ...
    'NormalizePathGains',true, ...
    'MaximumDopplerShift',30, ...
    'DopplerSpectrum',{doppler('Gaussian',0.6)}, ...
    'RandomStream','mt19937ar with seed', ...
    'Seed',22, ...
    'PathGainsOutputPort',true);


% Channel SNR
capacity=[];

channel = ones(symbol_count,1);


for step = power:1:lim
    
    a = 1;
    
    order = 2 ^ step;

    modem_name = sprintf('%1gQAM',order);

    index = 0 : order-1 ;

    qam_symbols = qammod(index,order);

    modulation = (1 / sqrt(sum(abs(qam_symbols).^2)/order) ) .* qam_symbols;

    for snr=SNR_dB
        
        % Generate some noise
        N0 = 1/(10^(snr/10));
        % Generate AWGN channel
        awgn_channel = comm.AWGNChannel('NoiseMethod',"Signal to noise ratio (Es/No)",'EsNo',snr);

        % Generate some random symbols
        symbols = ceil(length(modulation)*rand(1,symbol_count));

        % Generate the transmitted signal
        tx = modulation(symbols);
        
        % Generate the rayleigh channel signal
        %[tx, channel] = rayleigh_channel(tx');
        channel = complex(ones(symbol_count,1));
        rx = awgn_channel(tx);

        % Calculate the symbol probabilities
        probabilities0 = max( exp(-(abs( ones(length(modulation),1)*rx - modulation.'*channel').^2)/N0),  realmin );

        % Normalise the symbol probabilities
        probabilities = probabilities0 ./ (ones(length(modulation),1)*sum(probabilities0));

        % Calculate the capacity
        channel_capacity = log2(length(modulation)) + mean( sum(probabilities.*log2(probabilities)) );

        % Display the capacity
        %disp(['The channel capacity is ', num2str(channel_capacity), ' bits per channel use']);

        capacity(step,a)=channel_capacity;

        a=a+1;
    end
end

snr_dB = SNR_dB;
snr = 10.^(snr_dB/10);
C = log2(1+snr);

%capacity = load('-mat','capacity.mat').capacity;

figure
% Rayleigh fading channel CCMC capacity 

Legend = cell(9,1);

plot(SNR_dB,C);
Legend{1} = sprintf('CCMC Capacity');
grid on
xlabel('SNR (dB)');
ylabel('Channel Capacity (bit/s/Hz)');
title('The channel capacity of different modulation scheme over an AWGN and Rayleigh fading channel');
axis([SNR_dB(1) SNR_dB(end) 0 7])
hold on

plot(SNR_dB,capacity(1,:));
Legend{2}=sprintf('DCMC Capacity BPSK');
plot(SNR_dB,capacity(2,:));
Legend{3}=sprintf('DCMC Capacity QPSK');
plot(SNR_dB,capacity(3,:));
Legend{4}=sprintf('DCMC Capacity 8-QAM');
for i = 5:1:9
    plot(SNR_dB,capacity(i-1,:));
    Legend{i}=sprintf('DCMC Capacity %1gQAM',2^(i-1));
    hold on
end
legend(Legend,'Location','northwest');
grid on
xlabel('SNR (dB)');
ylabel('Channel Capacity (bit/s/Hz)');
title('Channel capacity of different modulation schemes over an AWGN channel');
axis([SNR_dB(1) SNR_dB(end) 0 7])

ylim([0 9]);

hold off

   