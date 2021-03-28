function [ leftSMRCSPFirst, leftSMRCSPLast, rightSMRCSPFirst, rightSMRCSPLast ] = mySMRCalculationCSP( filterB, filterA, EEG, leftEpochStartTime, rightEpochStartTime, trialTimeIdx, baselineIdx)

% filtering the signal 
for i=1:size(EEG, 1)
    eegFiltered(i, :) = filtfilt(filterB, filterA, EEG(i, :));
end

% extract epochs of left and right hand movements (including all channels)

% left
for epochID=1:length(leftEpochStartTime)
    indices = leftEpochStartTime(epochID) + trialTimeIdx;
    leftEpoches(:,:,epochID) = eegFiltered(:, indices);
end

% right
for epochID=1:length(rightEpochStartTime)
    indices = rightEpochStartTime(epochID) + trialTimeIdx;
    rightEpoches(:,:,epochID) = eegFiltered(:, indices);
end

% calculating the covarinace matrix for both classes (note CSP works on the
% raw signal, not the power).
% left class
SigmaL = zeros(size(leftEpoches, 1), size(leftEpoches, 1));
for epochID=1:length(leftEpochStartTime)
    SigmaL = SigmaL + cov(leftEpoches(:, :, epochID)');
end
SigmaL = SigmaL ./ length(leftEpochStartTime);

% right class
SigmaR = zeros(size(rightEpoches, 1), size(rightEpoches, 1));
for epochID=1:length(rightEpochStartTime)
    SigmaR = SigmaR + cov(rightEpoches(:, :, epochID)');
end
SigmaR = SigmaR ./ length(rightEpochStartTime);

% solving the generalized eigenvalue problem
% Note here SigmaL is the first input argumen. Here for simplicity, we take 
% the first column of W as the corresponding component for Left hand movement;
% and the last column of W is the corresponding component for right hand movement

[W,~] = eig(SigmaL, SigmaL + SigmaR);

% left class after CSP
for epochID = 1:size(leftEpoches,3)
    leftEpochesCSP(:,:,epochID) = W' * leftEpoches(:,:,epochID);
end
% right class after CSP
for epochID = 1:size(rightEpoches,3)
    rightEpochesCSP(:,:,epochID) = W' * rightEpoches(:,:,epochID);
end

% As noted above, we only take the first and the last component for left
% hand and right hand, respectively.
leftPowerSig = mean(leftEpochesCSP.^2, 3);
leftEpochesAvePowerCSPFirst = leftPowerSig(1, :);
leftEpochesAvePowerCSPLast = leftPowerSig(size(EEG, 1), :);

rightPowerSig = mean(rightEpochesCSP.^2, 3);
rightEpochesAvePowerCSPFirst = rightPowerSig(1, :);
rightEpochesAvePowerCSPLast = rightPowerSig(size(EEG, 1), :);

% find the power of the baseline range
leftBaselinePowerCSPFirst = mean(leftEpochesAvePowerCSPFirst(baselineIdx));
leftBaselinePowerCSPLast = mean(leftEpochesAvePowerCSPLast(baselineIdx));

rightBaselinePowerCSPFirst = mean(rightEpochesAvePowerCSPFirst(baselineIdx));
rightBaselinePowerCSPLast = mean(rightEpochesAvePowerCSPLast(baselineIdx));

% calculate the SMR of the trials (substract the baseline power, then
% normalized w.r.t. the baseline power)
leftSMRCSPFirst = ( (leftEpochesAvePowerCSPFirst - leftBaselinePowerCSPFirst) ./ leftBaselinePowerCSPFirst ) .* 100;
leftSMRCSPLast = ( (leftEpochesAvePowerCSPLast - leftBaselinePowerCSPLast) ./ leftBaselinePowerCSPLast ) .* 100;

rightSMRCSPFirst = ( (rightEpochesAvePowerCSPFirst - rightBaselinePowerCSPFirst) ./ rightBaselinePowerCSPFirst ) .* 100;
rightSMRCSPLast = ( (rightEpochesAvePowerCSPLast - rightBaselinePowerCSPLast) ./ rightBaselinePowerCSPLast ) .* 100;

end

