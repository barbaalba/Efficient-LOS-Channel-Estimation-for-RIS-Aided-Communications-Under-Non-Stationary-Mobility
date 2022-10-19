%% Initialization
numUE = 2;
RWL = 50; 
plt = true;
RIS_coor = [-25,0,2];
z_t = repelem(1.5,numUE,RWL+1);

%% Extract Azimuth-Elevation-Distance
% Random walk
[x_t,y_t] = randomwalk(numUE,RWL);

d_x = x_t-RIS_coor(1);
d_y = y_t-RIS_coor(2);
d_z = z_t-RIS_coor(3);
% Distance to RIS
d_t = sqrt(d_x.^2+d_y.^2+d_z.^2);
% Azimuth 
azimuth = d_y./d_x;
azimuth = atan(azimuth);

plotTrajectory(x_t,y_t,azimuth);

% Enlarge figure to full screen.
set(gcf, 'Units', 'Normalized', 'Outerposition', [0, 0.05, 1, 0.95]);
axis square;

