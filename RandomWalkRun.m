clear;clc;close all;
rng(13);
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
Speed = 0.5; % meter per second movement
plt = true; % To plot the trajectory
pltconf = 'discrete'; % 'continous' or 'discrete'
RIS_coor = [-Xmax,0,2]; % the RIS is located at the end of the room 
z_t = repelem(1.5,numUE,RWL+1); % Z-coordinate of the user(does not change)

%% Extract Azimuth-Elevation-Distance
% Random walk
[x_t,y_t] = randomwalk(numUE,RWL,Xmax,Ymax,Speed);
[azimuth,elevation,Cph,d_t] = ChanParGen(x_t,y_t,z_t,RIS_coor,lambda);
% report the near-field percentage
disp(['How much percentage in near field: ' num2str(sum(d_t < d_fraun,'all')/RWL/numUE*100)]);
if plt
    plotTrajectory(x_t,y_t,azimuth,elevation,pltconf,Xmax,Ymax,RIS_coor);
end
% uncomment it to save the required data
%save('fast20.mat','x_t','y_t','azimuth','elevation','Cph'); 

