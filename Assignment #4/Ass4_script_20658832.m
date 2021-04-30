% Assignment 4 script, SYDE544

% Winter 2021

% Student Name: Manthan Shah
% Student ID: 20658832

clear all;
close all;

rng('default'); % For reproducibility.

load Ass4_data;

% simulation duration
duration = 10;

% Problem 1a begins:
% generate the recruitment thresholds (RTE)

RR = 30;
a = log(RR) / 120;
for i = 1:size(MUAPs, 2)
   RTE(i) = exp(a * i);
end

% normalize, so the maximal RTE is 1
RTE = RTE / 30;

% Problem 1a ends


% Problem 1b begins:
% calculate the firing rate gain 
% minimal firing rate is 8 Hz, maximal rate is 35 Hz
for i = 1:size(MUAPs, 2)
    g(i) = (35 - 8) / (1.2 - RTE(i));
end

% get the firing rate at a certain excitation level
% Eqn (2) of Yao et al.
% Problem 1b ends


% Problem 2 begins
% Generate the firing timings for all units at the 9 neural drive levels

neuralDriveLvl = [3, 6, 9, 12, 15, 18, 21, 24, 27];

% Compute the firing rates for each MU at each neural drive level.
for i = 1:length(neuralDriveLvl)
    for j = 1:size(MUAPs, 2)
        
        % If the MU does not meet RTE then set 0 firing rate, as muscle
        % fibre will not get activated.
        if ( (neuralDriveLvl(i) / RR) - RTE(j) ) < 0
            firingRates(i, j) = 0;
        else
            firingRates(i, j) = (g(j) * ((neuralDriveLvl(i) / RR) - RTE(j))) + 8;
        end
    end
end

% Compute the firing timestamps for each MU at each neural drive level.
for i = 1:length(neuralDriveLvl)
    for j = 1:size(MUAPs, 2)
        t = 0;
        k = 1; % Variable to keep track of timestamp.
        
        % If the firing rate is 0, then set all timestamps to 0 for the MU
        % at current neural drive level.
        if firingRates(i, j) == 0
            firingsT(j, :, i) = 0;
        else
           meanFiringRate = 1 / firingRates(i, j);
           
           % Generate timestamps for the MU at current neural drive level
           % up to the simulation duration. Each ISI (interspike interval)
           % follows a Gaussian distribution.
           while t < duration
               if t == 0
                   t = normrnd(meanFiringRate, meanFiringRate * 0.15);
               else
                   t = firingsT(j, k - 1, i) + normrnd(meanFiringRate, meanFiringRate * 0.15);
               end
               
               % If timestamp exceeds simulation duration, drop it.
               if t < 10
                   firingsT(j, k, i) = t;
               end
               
               k = k + 1;
           end
        end
    end
end

% Problem 2 ends

% Problem 3 begins

% Sampling rate in Hz.
sr = 4096;

% Create impulse signal for each MU at each neural drive level.
impulseTrain = zeros(120, sr * duration, 9);

% Set impulses at each firing timestamp. This involves iterating each
% timestamp and converting it into an index in impulse train matrix.
for i = 1:length(neuralDriveLvl)
    for j = 1:size(MUAPs, 2)
        for k = 1:size(firingsT, 2)
            if firingsT(j, k, i) == 0
                break
            end
            
            % Translate firing timestamp into an index.
            idx = round(4096 * firingsT(j, k, i));
            
            % Set impulse at firing timestamp index.
            impulseTrain(j, idx, i) = 1;
        end
    end
end

% Create the MUAP trains.
for i = 1:length(neuralDriveLvl)
    for j = 1:size(MUAPs, 2)
        
        % MUAP train is a convolution operation between MUAP associated
        % with the MU, and corresponding firing timestamps.
        % X(t) = p(t) * u(t)
        MUAPtrain(j, :, i) = conv2(impulseTrain(j, :, i), MUAPs(:, j)');
    end
end

% simulated EMG signals as the summation of all MUAP trains
EMG = sum(MUAPtrain, 1);

% Problem 3 ends

% Problem 4 begins

% Calculate the SNR of EMG between 3s and 5s.
for i = 1:length(neuralDriveLvl)
    startIdx = 3 * 4096;
    endIdx = 5 * 4096;

    mu = mean(EMG(:, startIdx:endIdx, i).^2);
    v = var(EMG(:, startIdx:endIdx, i).^2);
    
    SNR_EMG(i) = mu^2 / v;
end

% Problem 4 ends


% Problem 5 begins    

% Creating force twitches.
for i = 1:length(neuralDriveLvl)
    for j = 1:size(MUAPs, 2)
        
        % Force twitch is a convolution operation between F associated
        % with the MU, and corresponding firing timestamps.
        forceTwt(j, :, i) = conv2(impulseTrain(j, :, i), F(:, j)');
    end
end

Force = sum(forceTwt, 1);

% Problem 5 ends            

% Problem 6 begins

% Calculate the SNR of Force between 3s and 5s.
for i = 1:length(neuralDriveLvl)
    startIdx = 3 * 4096;
    endIdx = 5 * 4096;
    
    mu = mean(Force(:, startIdx:endIdx, i));
    v = var(Force(:, startIdx:endIdx, i));
    
    SNR_Force(i) = mu^2 / v;
end

% Problem 6 ends

% Problem 7 begins

subplot(2, 1, 1);
plot(neuralDriveLvl, SNR_EMG);
title('SNR EMG');
xlabel('Neural Drive');
ylabel('SNR');

subplot(2, 1, 2);
plot(neuralDriveLvl, SNR_Force);
title('SNR Force');
xlabel('Neural Drive');
ylabel('SNR');

% problem 7 ends