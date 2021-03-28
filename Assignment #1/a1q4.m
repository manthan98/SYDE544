clear all;
close all;

vlims = [-3, 3];
wlims = [-3, 3];
v = linspace(vlims(1), vlims(2), 90);
w = linspace(wlims(1), wlims(2), 90)';

epsilon = 0.1;
b0 = 2;
b1 = 1.5;
Iapp = 0;

dummyV = repmat(v, length(w), 1);
dummyW = repmat(w, 1, length(v));

dv = dummyV - dummyV.^3/3 - dummyW + Iapp;
dw = epsilon * (b0 + b1.*dummyV - dummyW);

% v-nullcline:
vNullcline = v - v.^3/3 + Iapp;

% w-nullcine:
wNullcline = b0 + b1.*v;

% Draw the vector field
startPoint = 1;
interval = 3;
figure;
qObj(1) = quiver(dummyV(startPoint:interval:end,startPoint:interval:end), dummyW(startPoint:interval:end,startPoint:interval:end),...
    dv(startPoint:interval:end,startPoint:interval:end), dw(startPoint:interval:end,startPoint:interval:end));

% Plotting the nullclines
vLine(1) = line(v, vNullcline, 'color', 'g', 'linewidth', 2);
wLine(2) = line(w, wNullcline, 'color', 'r', 'linewidth', 2);
set(gca,'xlim', vlims, 'ylim', wlims);
xlabel('V');
ylabel('W');

qObj.AutoScaleFactor = 0.9;

% Start the simulation
points = 1;
while points
    startPos = ginput(1);
    if exist('hNullclineTraj','var')
        delete(hNullclineTraj);
        delete(hCurrentPoint);
    end
    
    maxStep = 5e4;
    minMove = 1e-6;
    vw = startPos;
    
    hNullclineTraj = line(vw(1), vw(2), 'color', 'b', 'linewidth', 1.5);
    hCurrentPoint = line(vw(1), vw(2), 'color', 'r', 'linewidth', 1, 'marker', 'o', 'markerfacecolor', 'red', 'markersize', 4);
    
    dt = 0.01;
    for nStep = 1:1:maxStep
        c_dv = (vw(nStep, 1) - vw(nStep, 1)^3/3 - vw(nStep, 2) + Iapp) * dt;
        c_dw =  epsilon * (b0 + b1*vw(nStep, 1) - vw(nStep, 2)) * dt;
        
        vw(nStep+1,1) = vw(nStep,1) + c_dv;
        vw(nStep+1,2) = vw(nStep,2) + c_dw;
        set(hNullclineTraj, 'xdata', vw(:,1), 'ydata', vw(:,2));
        set(hCurrentPoint, 'xdata', vw(end,1), 'ydata', vw(end,2));
        
        pause(1e-4);
        if sqrt((vw(nStep+1,1)-vw(nStep,1))^2 + (vw(nStep+1,2) - vw(nStep,2))^2) < minMove
            break;
        end
    end
    points = points + 1;
end




