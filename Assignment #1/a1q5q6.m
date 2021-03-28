clear all;
close all;

isQ6 = false; % Set this to 'true' for question #6

vlims = [-3, 3];
wlims = [-3, 3];
v = linspace(vlims(1), vlims(2), 150);
w = linspace(wlims(1), wlims(2), 150)';

epsilon = 0.1;
b0 = 2;
b1 = 1.5;
Iapp = 0;
if isQ6
    Iapp = 1;
end

dummyV = repmat(v, length(w), 1);
dummyW = repmat(w, 1,length(v));

dx = dummyV - dummyV.^3/3 - dummyW + Iapp;
dy = epsilon * (b0 + b1.*dummyV - dummyW);

% v-nullcline:
vNullcline = v - v.^3/3 + Iapp;

% w-nullcine:
wNullcline = b0 + b1.*v;

% Draw the vector field
startPoint = 1;
interval = 5;
figure,subplot(1,2,1);
qObj(1) = quiver(dummyV(startPoint:interval:end,startPoint:interval:end),dummyW(startPoint:interval:end,startPoint:interval:end),...
    dx(startPoint:interval:end,startPoint:interval:end),dy(startPoint:interval:end,startPoint:interval:end));
set(qObj(1),'ShowArrowHead','off');

% Plotting the nullclines
vLine(1) = line(v,vNullcline,'color','g','linewidth',2);
wLine = line(w,wNullcline,'color','r','linewidth',2);
axs1 = get(gcf,'children');
xlabel('V');
ylabel('W');
set(axs1,'xlim',vlims,'ylim',wlims,'fontsize',15);
qObj.AutoScaleFactor = 0.9;

positions = [-2 -2/3; -1 -2/3; 0 -2/3];

% Change the position index for different simulation
startPos = positions(1, :);

% Start the simulation
if exist('hUTraj','var')
    delete(hNullclineTraj);
    delete(hCurrentPoint);
    delete(hUTraj);
end

maxStep = 5e4;
minMove = 1e-6;
vw = [startPos(1), startPos(2)];
hNullclineTraj = line(vw(1),vw(2),'color','b','linewidth',1.5);
hCurrentPoint = line(vw(1),vw(2),'color','r','linewidth',1,'marker','o','markerfacecolor','red','markersize',4);
vAxes = subplot(1,2,2);
xlabel('t');
ylabel('V')
set(vAxes,'fontsize',15);
hUTraj = line(0,vw(1),'color','b','linewidth',1);

dt = 0.1;
for nStep = 1:maxStep
    c_dv = (vw(nStep, 1) - vw(nStep, 1)^3/3 - vw(nStep, 2) + Iapp) * dt;
    c_dw = epsilon * (b0 + b1*vw(nStep, 1) - vw(nStep, 2)) * dt;
    
    vw(nStep+1,1) = vw(nStep,1) + c_dv;
    vw(nStep+1,2) = vw(nStep,2) + c_dw;
    if rem(nStep, 10) == 0
        set(hNullclineTraj,'xdata',vw(:,1),'ydata',vw(:,2));
        set(hCurrentPoint,'xdata',vw(end,1),'ydata',vw(end,2));
        set(hUTraj,'ydata',vw(:,1),'xdata',(0:nStep)/dt);
        pause(1e-5);
    end
    if sqrt((vw(nStep+1,1)-vw(nStep,1))^2 + (vw(nStep+1,2) - vw(nStep,2))^2) < minMove
        break;
    end
end
