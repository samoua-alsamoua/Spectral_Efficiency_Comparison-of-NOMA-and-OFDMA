% MATLAB Code for Spectral Efficiency Comparison: NOMA vs OFDMA

clc;
clear;
close all;

% Parameters
numUsers = 4;               % Number of users
numSubcarriers = 64;         % Number of subcarriers (must be divisible by numUsers)
modulationOrder = [4, 16];  % Modulation orders (QPSK, 16-QAM)
bandwidth = 10e6;           % Bandwidth in Hz (10 MHz)
snr_dB = 0:2:20;            % SNR range in dB
numIter = 100;              % Number of iterations for averaging

% Ensure numSubcarriers is divisible by numUsers
if mod(numSubcarriers, numUsers) ~= 0
    error('Number of subcarriers must be divisible by the number of users.');
end

% Initialize metrics
bitsPerSymbol_NOMA = zeros(length(snr_dB), 1);
bitsPerSymbol_OFDMA = zeros(length(snr_dB), 1);
bitsPerHz_NOMA = zeros(length(snr_dB), 1);
bitsPerHz_OFDMA = zeros(length(snr_dB), 1);
oobEmissions_NOMA = zeros(length(snr_dB), 1);
oobEmissions_OFDMA = zeros(length(snr_dB), 1);

% Main loop for SNR
for snrIdx = 1:length(snr_dB)
    snr = snr_dB(snrIdx);
    bps_NOMA = 0;
    bps_OFDMA = 0;
    bph_NOMA = 0;
    bph_OFDMA = 0;
    oob_NOMA = 0;
    oob_OFDMA = 0;
    
    for iter = 1:numIter
        % Generate random data for each user
        data_NOMA = randi([0, modulationOrder(2) - 1], numSubcarriers, numUsers);
        data_OFDMA = randi([0, modulationOrder(2) - 1], numSubcarriers, numUsers);
        
        % Modulate data using QPSK and 16-QAM
        modData_NOMA = qammod(data_NOMA(:), modulationOrder(2), 'UnitAveragePower', true);
        modData_OFDMA = qammod(data_OFDMA(:), modulationOrder(2), 'UnitAveragePower', true);
        
        % Reshape modulated data
        modData_NOMA = reshape(modData_NOMA, numSubcarriers, numUsers);
        modData_OFDMA = reshape(modData_OFDMA, numSubcarriers, numUsers);
        
        % NOMA: Superposition coding
        powerAllocation = linspace(1, 0.5, numUsers); % Power allocation for users
        nomaSignal = sum(modData_NOMA .* powerAllocation, 2);
        
        % OFDMA: Assign subcarriers to users
        ofdmaSignal = zeros(numSubcarriers, 1);
        subcarriersPerUser = numSubcarriers / numUsers; % Subcarriers per user
        for u = 1:numUsers
            startIdx = (u-1) * subcarriersPerUser + 1;
            endIdx = u * subcarriersPerUser;
            ofdmaSignal(startIdx:endIdx) = modData_OFDMA(startIdx:endIdx, u);
        end
        
        % Add AWGN noise
        nomaRx = awgn(nomaSignal, snr, 'measured');
        ofdmaRx = awgn(ofdmaSignal, snr, 'measured');
        
        % NOMA Receiver: SIC (Successive Interference Cancellation)
        nomaData = qamdemod(nomaRx, modulationOrder(2), 'UnitAveragePower', true);
        
        % OFDMA Receiver: Extract user data
        ofdmaData = qamdemod(ofdmaRx, modulationOrder(2), 'UnitAveragePower', true);
        
        % Calculate Bits per Symbol
        bps_NOMA = bps_NOMA + mean(log2(modulationOrder(2)) * numUsers);
        bps_OFDMA = bps_OFDMA + mean(log2(modulationOrder(2)) * numUsers);
        
        % Calculate Bits per Hz
        bph_NOMA = bph_NOMA + (mean(log2(modulationOrder(2)) * numUsers) / (bandwidth / numSubcarriers));
        bph_OFDMA = bph_OFDMA + (mean(log2(modulationOrder(2)) * numUsers) / (bandwidth / numSubcarriers));
        
        % Calculate Out-of-Band Emissions (simplified as power spectral density)
        [psdNOMA, ~] = pwelch(nomaSignal, [], [], [], bandwidth);
        [psdOFDMA, ~] = pwelch(ofdmaSignal, [], [], [], bandwidth);
        oob_NOMA = oob_NOMA + mean(psdNOMA(1:10)); % Average OOB emissions
        oob_OFDMA = oob_OFDMA + mean(psdOFDMA(1:10)); % Average OOB emissions
    end
    
    % Average metrics
    bitsPerSymbol_NOMA(snrIdx) = bps_NOMA / numIter;
    bitsPerSymbol_OFDMA(snrIdx) = bps_OFDMA / numIter;
    bitsPerHz_NOMA(snrIdx) = bph_NOMA / numIter;
    bitsPerHz_OFDMA(snrIdx) = bph_OFDMA / numIter;
    oobEmissions_NOMA(snrIdx) = oob_NOMA / numIter;
    oobEmissions_OFDMA(snrIdx) = oob_OFDMA / numIter;
end

% Plot Results
figure;
subplot(3, 1, 1);
plot(snr_dB, bitsPerSymbol_NOMA, 'bo-', 'LineWidth', 2);
hold on;
plot(snr_dB, bitsPerSymbol_OFDMA, 'rs-', 'LineWidth', 2);
grid on;
xlabel('SNR (dB)');
ylabel('Bits per Symbol');
legend('NOMA', 'OFDMA');
title('Bits per Symbol Comparison');

subplot(3, 1, 2);
plot(snr_dB, bitsPerHz_NOMA, 'bo-', 'LineWidth', 2);
hold on;
plot(snr_dB, bitsPerHz_OFDMA, 'rs-', 'LineWidth', 2);
grid on;
xlabel('SNR (dB)');
ylabel('Bits per Hz');
legend('NOMA', 'OFDMA');
title('Spectral Efficiency (Bits per Hz)');

subplot(3, 1, 3);
plot(snr_dB, oobEmissions_NOMA, 'bo-', 'LineWidth', 2);
hold on;
plot(snr_dB, oobEmissions_OFDMA, 'rs-', 'LineWidth', 2);
grid on;
xlabel('SNR (dB)');
ylabel('Out-of-Band Emissions');
legend('NOMA', 'OFDMA');
title('Out-of-Band Emissions Comparison');