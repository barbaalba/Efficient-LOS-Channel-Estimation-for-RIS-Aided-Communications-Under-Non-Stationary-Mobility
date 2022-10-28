clear;clc;close all;
%% Scanario Initialization
freq = 28e9; % Central frequency
lambda = physconst('LightSpeed') / freq; % Wavelength
M_H = 64; M_V = 64; M = M_H * M_V;
RISspacing = 1/2;
disp(['X-axis size is ' num2str(M_H*RISspacing*lambda) ' (m)']);
disp(['Y-axis size is ' num2str(M_V*RISspacing*lambda) ' (m)']);
d_fraun = 4*(M_H*RISspacing*lambda)^2/lambda/10;
disp(['Fraunhofer distance is ' num2str(d_fraun) ' (m)']);
d_BSRIS = 30;
Xmax = 5; % Room size = [-Xmax, Xmax]
Ymax = 5; % Room size = [-Ymax, Ymax]
numUE = 1; % number of users
RWL = 200;  % Length of the random walk
Speed = 0.25; % meter per second movement
plt = true; % To plot the trajectory
pltconf = 'discrete'; % 'continous' or 'discrete'
RIS_coor = [-Xmax,0,2]; % the RIS is located at the end of the room 
z_t = repelem(1.5,numUE,RWL+1); % Z-coordinate of the user(does not change)
% channel parameters (LOS mmWave) alpha + 10*beta*log10(d)
alpha = 61.4; % path loss reference in 1 (m) distance
beta = 1.46;
plEval = @(d) alpha + 10*beta*log10(d);
noisepow = -96; % [dBm]
txpow = 0; % [dBm]
%% Extract Azimuth-Elevation-Distance
% Random walk
[x_t,y_t] = randomwalk(numUE,RWL,Xmax,Ymax,Speed);

d_x = x_t-RIS_coor(1);
d_y = y_t-RIS_coor(2);
d_z = z_t-RIS_coor(3);
% Distance to RIS
d_t = sqrt(d_x.^2+d_y.^2+d_z.^2);
Cph = exp(-1i*2*pi*d_t/lambda);
% report the near-field percentage
sum(d_t < d_fraun,'all')/RWL/numUE*100
% Azimuth 
azimuth = d_y./d_x;
azimuth = atan(azimuth);
elevation =  d_z ./ d_t;
elevation = asin(elevation);
if plt
    plotTrajectory(x_t,y_t,azimuth,elevation,pltconf,Xmax,Ymax,RIS_coor);
end
% uncomment it to save the required data
save('slow.mat','x_t','y_t','azimuth','elevation','Cph'); 
%% Direct Channel model 
g_d = -130; % the NLOS direct channel gain
h_d = sqrt(db2pow(g_d)) * (randn+1i*randn); % direct channel complex gain
%% End-End channel 
% channel between RIS and BS
PLdB = plEval(d_BSRIS);
RIS_AOA = pi/6; % Azimuth of arrival at RIS from BS
RIS_EOA = 0; % Elevation of arrival to RIS from BS (They are in the same level)
g = sqrt(db2pow(PLdB))*exp(1i*2*pi/lambda*d_BSRIS)*...
    UPA_Evaluate(lambda,M_V,M_H,RIS_AOA,RIS_EOA,RISspacing,RISspacing);
% channel between RIS and UE and the algorithm
UEPL = plEval(d_t);
disp(['End to End path loss is ' num2str(mean(UEPL,'all')+PLdB) '[dB]']);
h = zeros(M_H*M_V,RWL+1,numUE);
for k=1:numUE
    h(:,:,k) = sqrt(db2pow(-UEPL(k,:))) .*UPA_Evaluate(lambda,M_V,M_H,azimuth(k,:),elevation(k,:),RISspacing,RISspacing);
end

% Enlarge figure to full screen.
% set(gcf, 'Units', 'Normalized', 'Outerposition', [0, 0.05, 1, 0.95]);
% axis square;

