clear all;
close all;

% To grader: I have completed all questions, please feel free to give me
% feedback on some/all questions if possible. I had completed the
% entire assignment before the announcement to only do question 1. Thanks.

% Student Name: Manthan Shah
% Student ID: 20658832

% P-QRS-T duration and axis
% duration
P_duration = 110; % ms
PR_duration = 40; % ms
QRS_duration = 100; % ms
ST_duration = 90; %ms
T_duration = 150; % ms

% wave axis
PWave_axis = pi*55/180;
QRS_axis = 45*pi/180;
T_axis = pi*45/180;

sr = 1; % kHz

% generate QRS Loop
% Beginning of Question 1
% Question 1 with Eqn (1) and proper rotation
%
a = 4;
b = 1/3;
theta = linspace(0, 2*pi, sr * QRS_duration);
QRS_Loop = [(a.*(cos(theta) - 1).*cos(theta))', (b.*a.*(cos(theta) - 1).*sin(theta))'];
QRS_Loop = QRS_Loop*([[cos(QRS_axis), -sin(QRS_axis)];[sin(QRS_axis), cos(QRS_axis)]]); % Rotate
% End of Question 1


% generate P Loop
% Beginning of Question 2
% Question 2 with Eqn (2) and proper rotation
%
a = 1;
b = 1/3;
theta = linspace(-pi, pi, sr * P_duration);
P_Loop = [(a + a.*cos(theta))', (b.*sin(theta))'];
P_Loop = P_Loop*([[cos(PWave_axis), -sin(PWave_axis)];[sin(PWave_axis), cos(PWave_axis)]]); % Rotate
% End of Question 2

% generate T Loop
% Beginning of Question 3
% Question 3 with Eqn (2) and proper rotation
%
a = 1.5;
b = 0.5;
theta = linspace(-pi, pi, sr * T_duration);
T_Loop = [(a + a.*cos(theta))', (b.*sin(theta))'];
T_Loop = T_Loop*([[cos(T_axis), -sin(T_axis)];[sin(T_axis), cos(T_axis)]]); % Rotate
% End of Question 3 


% Generate the signals in six limb lead (only I, II, III, aVL, aVR, and aVF)

% Beginning of question 4
% angles of leads (note definitions given in class are in the clock-wise direction)
Lead_Angles = [0, -pi/3, -2*pi/3, pi/6, 5*pi/6, -pi/2]; % [Phase I, Phase II, Phase III, aVL, aVR, aVF]

waveform_duration = (P_duration + PR_duration + QRS_duration + ST_duration + T_duration + 100) * sr;
num_repeats = 3;
ECG = zeros(waveform_duration * num_repeats, 6);

for ecgLeads = 1:6
    % https://www.mathworks.com/help/matlab/ref/pol2cart.html#bu5a0kc-theta
    [x, y] = pol2cart(Lead_Angles(ecgLeads), 1); % Vector in direction of lead angle.
    proj_vec = [x, y];
    
    idx = 1;
    
    for i = 1:num_repeats
        for j = 1:length(P_Loop)
            vec = P_Loop(j, :);
            
            % https://www.mathworks.com/matlabcentral/answers/2216-projecting-a-vector-to-another-vector
            projected = (dot(vec,proj_vec)/norm(proj_vec)^2)*proj_vec;
            projected = (projected(1)^2 + projected(2)^2)^0.5;

            % https://www.mathworks.com/matlabcentral/answers/101590-how-can-i-determine-the-angle-between-two-vectors-in-matlab
            CosTheta = max(min(dot(vec,proj_vec)/(norm(vec)*norm(proj_vec)),1),-1);
            ThetaInDegrees = real(acosd(CosTheta));

            if ThetaInDegrees > 90 || ThetaInDegrees < -90
                projected = -1 * projected;
            end

            ECG(idx, ecgLeads) = projected;
            idx = idx + 1;
        end
        
        idx = idx + PR_duration;
        
        for j = 1:length(QRS_Loop)
            vec = QRS_Loop(j, :);
            projected = (dot(vec,proj_vec)/norm(proj_vec)^2)*proj_vec;
            projected = (projected(1)^2 + projected(2)^2)^0.5;

            CosTheta = max(min(dot(vec,proj_vec)/(norm(vec)*norm(proj_vec)),1),-1);
            ThetaInDegrees = real(acosd(CosTheta));

            if ThetaInDegrees > 90 || ThetaInDegrees < -90
                projected = -1 * projected;
            end

            ECG(idx, ecgLeads) = projected;
            idx = idx + 1;
        end
        
        idx = idx + ST_duration - 1;
        
        for j = 1:length(T_Loop)
            vec = T_Loop(j, :);
            projected = (dot(vec,proj_vec)/norm(proj_vec)^2)*proj_vec;
            projected = (projected(1)^2 + projected(2)^2)^0.5;

            CosTheta = max(min(dot(vec,proj_vec)/(norm(vec)*norm(proj_vec)),1),-1);
            ThetaInDegrees = real(acosd(CosTheta));

            if ThetaInDegrees > 90 || ThetaInDegrees < -90
                projected = -1 * projected;
            end

            ECG(idx, ecgLeads) = projected;
            idx = idx + 1;
        end
        
        idx = idx + 100;
    end
end

% visualize the simulated signals
figure;

% Projection of loops onto Phase I lead at 0 degrees.
% P-wave: larger peak magnitude than all but phase II and aVF ECG leads.
%   - P-Loop is at 55 degree angle (CW) wrt phase I ECG lead. 
%   - The wave is negative initially as depolarization vector runs in
%   opposite direction of phase I lead, but positive for rest of the wave 
%   after.
%   - The magnitude is smaller than QRS and T waves since loop major and
%   minor axes are smaller in magnitude.
% QRS-wave: Q wave peak magnitude is smaller than all but aVL ECG lead, R 
% wave smaller than phase II ECG lead, similar to aVF, and larger than all 
% other leads, and S wave is larger than all but aVL and aVR ECG leads.
%   - QRS-loop is at 45 degree angle (CW) wrt phase I. 
%   - Q wave dips negative since depolarization vector runs in opposite 
%   direction to phase I, R wave goes positive as depolarization vector 
%   turns in direction of phase I, and then S wave dips negative as 
%   depolarization vector runs in opposite direction of phase I lead again.
% T-wave: larger than all other leads other than phase II ECG lead.
%   - T-loop is at 45 degree angle (CW) wrt phase I. 
%   - T wave dips negative briefly, as depolarization vector runs in 
%   opposite direction of phase I, but then it goes positive as 
%   depolarization vector turns in direction of phase I.
subplot(3,2,1), plot(ECG(:,1)), title('Phase I ECG'), xlabel('Time (ms)'), ylabel('mV');

% Projection of loops onto Phase II lead at 60 degrees (CW).
%   - Large magnitudes of all peaks due to close orientation of loops to
%   phase II lead.
% P-wave: largest peak magnitude than all other leads.
%   - P-Loop is at 5 degree angle (CCW) wrt phase II. 
%   - It is always positive, thus the depolarization vector direction is 
%   always in same direction as phase II lead.
% QRS-wave: Q wave is more negative than phase III, aVR, aVF but less than 
% phase I, aVL ECG leads. R wave is larger than all other ECG leads. S wave
% is smaller than all ECG leads except phase III and aVF.
%   - QRS-Loop is at 15 degree angle (CCW) wrt phase II lead.
%   - Q wave dips negative since depolarization vector runs in opposite 
%   direction to phase I, R wave goes positive as depolarization vector 
%   turns in direction of phase I, and then S wave dips negative as 
%   depolarization vector runs in opposite direction of phase I lead again.
% T-wave: larger peak magnitude than all other leads.
%   - T-Loop is at 15 degree angle (CCW) wrt phase II lead. 
%   - Wave is positive as depolarization vector runs in same direction as 
%   phase II lead at all times.
subplot(3,2,3),plot(ECG(:,2)), title('Phase II ECG'), xlabel('Time (ms)'), ylabel('mV');

% Projection of loops onto Phase III lead at 120 degrees (CW).
%  - Smaller magnitudes (compared to other leads) as difference between 
%   direction of loops and phase III lead is large.
% P-wave: smaller peak magnitude than all but aVL and aVR ECG leads.
%   - P-Loop is at 65 degree angle (CCW) wrt phase III. 
%   - It is positive for the majority of the peak as the depolarization 
%   vector is in direction of phase III lead, and then dips slightly 
%   negative at the end due to opposite depolarization vector direction 
%   with phase III lead.
% QRS-wave: Q wave larger than all except aVR ECG lead, R wave is similar
% to aVL, smaller than all others except aVR, and S wave is similar to aVF
% and smaller than all other ECG leads.
%   - QRS-Loop is at 75 degree angle (CCW) wrt phase III.
%   - Q wave dips negative since depolarization vector runs in opposite 
%   direction to phase III, R wave goes positive as depolarization vector 
%   turns in direction of phase III, and then S wave dips negative as 
%   depolarization vector runs in opposite direction of phase III lead again.
% T-wave: peak magnitude similar to aVL ECG lead, and smaller than all 
% others except aVR ECG lead.
%   - T-Loop is at 75 degree angle (CCW) wrt phase III. 
%   - Initially, the wave is positive as depolarization vector is in same 
%   direction as phase III lead, and then the wave dips negative as 
%   depolarization vector runs in opposite direction of phase III lead.
subplot(3,2,5),plot(ECG(:,3)), title('Phase III ECG'), xlabel('Time (ms)'), ylabel('mV');

% Projection of loops onto aVL lead at 330 degrees (CW).
%   - Smaller magnitudes (compared to other leads) as difference between 
%   direction of loops and aVL lead is large.
% P-wave: smaller peak magnitude than all leads, except aVR.
%   - P-Loop is at 85 degree angle (CW) wrt aVL lead. 
%   - It dips negative initially, as depolarization vector runs in 
%   direction away from lead, and then the wave goes positive as the 
%   depolarization vector approaches in the same direction as lead.
% QRS-wave: Q wave peak magnitude is smaller than all other leads, R wave
% is similar to phase III and smaller than all others except aVR, and S
% wave is larger than all except aVR ECG lead.
%   - QRS-Loop is at 75 degree angle (CW) wrt aVL.
%   - Q wave dips negative since depolarization vector runs in opposite 
%   direction to aVL, R wave goes positive as depolarization vector 
%   turns in direction of aVL, and then S wave dips negative as 
%   depolarization vector runs in opposite direction of aVL lead again.
% T-wave: similar peak magnitude to phase III ECG lead but smaller than all
% except aVR lead.
%   - T-Loop is at 75 degree angle (CW) wrt aVL.
%   - T wave dips negative briefly, as depolarization vector runs in 
%   opposite direction of aVL, but then it goes positive as 
%   depolarization vector turns in direction of aVL.
subplot(3,2,2),plot(ECG(:,4)), title('aVL ECG'), xlabel('Time (ms)'), ylabel('mV');

% Projection of loops onto aVR lead at 210 degrees (CW).
%   - Large magnitudes in the negative direction (compared to other leads)
%   due to large difference between direction of loops and aVR lead.
%   However, the directions are oppositely parallel, leading to large
%   magnitudes.
% P-wave: smaller peak magnitude than all other leads.
%   - P-Loop is at 205 degree angle (CW) wrt aVR lead. 
%   - Wave is negative as depolarization vector runs in opposite direction 
%   of aVR lead at all times.
% QRS-wave: Q wave magnitude is larger than all other leads, R wave is
% smaller than all other leads, and S wave is larger than all other leads.
%   - QRS-Loop is at 195 degree angle (CW) wrt aVR lead. 
%   - Q wave is positive as the depolarization vector runs in same direction 
%   as aVR lead initially, then R wave is negative as the depolarization 
%   vector runs in opposite direction of aVR, and finally the S wave 
%   becomes positive as depolarization vector runs in same direction as 
%   aVR lead again.
% T-wave: smaller peak magnitude than all other leads.
%   - T-Loop is at 195 degree angle (CW) wrt aVR lead. 
%   - Wave is negative as depolarization vector runs in opposite direction 
%   of aVR lead at all times.
subplot(3,2,4),plot(ECG(:,5)), title('aVR ECG'), xlabel('Time (ms)'), ylabel('mV');

% Projection of loops onto aVF lead at 90 degrees (CW).
%   - Large magnitudes (compared to other leads) as difference between
%   direction of loops and aVF lead is small.
% P-wave: larger peak magnitude than all others except phase II ECG lead.
%   - P-Loop is at 35 degree angle (CW) wrt aVF lead.
%   - It is always positive, thus the depolarization vector direction is 
%   always in same direction as aVF lead.
% QRS-wave: Q wave has larger magnitude than all other leads except phase
% III and aVR, R wave is similar to phase I and larger than all others
% except phase II, and S wave is smaller than all ECG leads except phase
% III ECG lead.
%   - QRS-Loop is 45 degree angle (CCW) wrt aVF lead.
%   - Q wave dips negative since depolarization vector runs in opposite 
%   direction to aVF, R wave goes positive as depolarization vector 
%   turns in direction of aVF, and then S wave dips negative as 
%   depolarization vector runs in opposite direction of aVF lead again.
% T-wave: T-wave is is similar to phase I ECG, smaller than phase II ECG, 
% and larger than all other ECG leads.
%   - T-Loop is 45 degree angle (CW) wrt aVF lead.
%   - Wave is positive as depolarization vector runs in same direction as 
%   aVF lead at all times.
subplot(3,2,6),plot(ECG(:,6)), title('aVF ECG'), xlabel('Time (ms)'), ylabel('mV');

% End of Question 4