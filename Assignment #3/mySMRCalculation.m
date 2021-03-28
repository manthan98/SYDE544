function [ C3LeftSMR, C4LeftSMR, C3RightSMR, C4RightSMR ] = ...
    mySMRCalculation( filterB, filterA, C3, C4, leftEpochStartTime, rightEpochStartTime, trialTimeIdx, baselineIdx)


% filtering the signal 
C3Filtered = filtfilt(filterB, filterA, C3);
C4Filtered = filtfilt(filterB, filterA, C4);

% extract epochs of left and right hand movements
% power of the left epoches (how to get the power?)
% 
for epochID=1:length(leftEpochStartTime)
    indices = trialTimeIdx + leftEpochStartTime(epochID);
    C3LeftEpoches(epochID,:) = C3Filtered(indices).^2;
    C4LeftEpoches(epochID,:) = C4Filtered(indices).^2;
end

% power of the right epoches
for epochID=1:length(rightEpochStartTime)
    indices = trialTimeIdx + rightEpochStartTime(epochID);
    C3RightEpoches(epochID,:) = C3Filtered(indices).^2;
    C4RightEpoches(epochID,:) = C4Filtered(indices).^2;
end

% get the ave. power of all trials - mean of each column
C3LeftEpochesAvePower = mean(C3LeftEpoches);
C4LeftEpochesAvePower = mean(C4LeftEpoches);

C3RightEpochesAvePower = mean(C3RightEpoches);
C4RightEpochesAvePower = mean(C4RightEpoches);

% find the power of the baseline range
C3LeftBaselinePower = mean(C3LeftEpochesAvePower(baselineIdx));
C4LeftBaselinePower = mean(C4LeftEpochesAvePower(baselineIdx));

C3RightBaselinePower = mean(C3RightEpochesAvePower(baselineIdx));
C4RightBaselinePower = mean(C4RightEpochesAvePower(baselineIdx));

% calculate the SMR of the trials (substract the baseline, then normalized
% w.r.t. the baseline
C3LeftSMR = ((C3LeftEpochesAvePower - C3LeftBaselinePower) ./ C3LeftBaselinePower) .* 100;
C4LeftSMR = ((C4LeftEpochesAvePower - C4LeftBaselinePower) ./ C4LeftBaselinePower) .* 100;

C3RightSMR = ((C3RightEpochesAvePower - C3RightBaselinePower) ./ C3RightBaselinePower) .* 100;
C4RightSMR = ((C4RightEpochesAvePower - C4RightBaselinePower) ./ C4RightBaselinePower) .* 100;

end

