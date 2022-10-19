%% Initialization
Xmax = 25; % Room size = [-Xmax, Xmax]
Ymax = 25; % Room size = [-Ymax, Ymax]
numUE = 2; % number of users
RWL = 500;  % Length of the random walk
plt = true; % To plot the trajectory
pltconf = 'discrete'; % 'continous' or 'discrete'
RIS_coor = [-Ymax,0,2]; % the RIS is located at the end of the room 
z_t = repelem(1.5,numUE,RWL+1); % Z-coordinate of the user(does not change)

%% Extract Azimuth-Elevation-Distance
% Random walk
[x_t,y_t] = randomwalk(numUE,RWL,Xmax,Ymax);

d_x = x_t-RIS_coor(1);
d_y = y_t-RIS_coor(2);
d_z = z_t-RIS_coor(3);
% Distance to RIS
d_t = sqrt(d_x.^2+d_y.^2+d_z.^2);
% Azimuth 
azimuth = d_y./d_x;
azimuth = atan(azimuth);

plotTrajectory(x_t,y_t,azimuth,pltconf);

% Enlarge figure to full screen.
set(gcf, 'Units', 'Normalized', 'Outerposition', [0, 0.05, 1, 0.95]);
axis square;

