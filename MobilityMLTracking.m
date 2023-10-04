% The code is intended to understand the performance of the proposed
% tracking sheme under mobility scenario 
close all; clc; clear; 
freq = 28e9; % Central frequency
lambda = physconst('LightSpeed') / freq; % Wavelength
rng(4);
% Mobility Config (Room Size, User speed, RIS coordination)
Speed = 0.1; numUE = 1; RWL = 1000;
Xmax = 5; Ymax = 5; z_t = repelem(1.5,numUE,RWL+1); randchan = true;
RIS_coor = [-Xmax,0,2];
% Channel parameters
if ~randchan 
    load("fast.mat","azimuth","Cph","elevation");
else
    [x_t,y_t] = randomwalk(numUE,RWL,Xmax,Ymax,Speed);
    [azimuth,elevation,Cph,d_t] = ChanParGen(x_t,y_t,z_t,RIS_coor,lambda);
end

% Tracking config
Prep = 50; % period of channel estimation

%UPA Element configuration
M_H = 8; M_V = 8; M = M_H*M_V;
elementspacing = 1/4; %In wavelengths

%Set the SNR
SNRdB_pilot = 0;
SNR_pilot = db2pow(SNRdB_pilot);

SNRdB_data = -10;
SNR_data = db2pow(SNRdB_data);

%Select angle to the base station (known value)
varphi_BS = -pi/6;
theta_BS = 0;

%Generate channels
h = UPA_Evaluate(lambda,M_V,M_H,varphi_BS,theta_BS,elementspacing,elementspacing);
Dh = diag(h);
Dh_angles = diag(h./abs(h)); 

nbrOfAngleRealizations = RWL;
nbrOfNoiseRealizations = 100;


%Save the rates achieved at different iterations of the algorithm
capacity = zeros(1,nbrOfAngleRealizations);
rate_proposed = zeros(M-1,nbrOfAngleRealizations,nbrOfNoiseRealizations);
rate_LS = zeros(M-1,nbrOfAngleRealizations,nbrOfNoiseRealizations);


%Create a uniform grid of beams (like a DFT matrix) to be used at RIS
ElAngles = asin((-M_V/2:1:M_V/2-1)*2/M_V);
AzAngles = asin((-M_H/2:1:M_H/2-1)*2/M_H);
beamAngles = zeros(M,2); % Elevation-Azimuth pair
beamresponses = zeros(M,M);
for i = 1:length(ElAngles)
    beamAngles ((i-1)*length(AzAngles)+1: i*length(AzAngles),:) = [repelem(ElAngles(i),length(AzAngles),1) AzAngles'];
    beamresponses(:,(i-1)*length(AzAngles)+1: i*length(AzAngles)) = UPA_Evaluate(lambda,M_V,M_H,AzAngles,repelem(ElAngles(i),1,length(AzAngles)),elementspacing,elementspacing);
end

% Define a fine grid of angle directions to analyze when searching for angle of arrival
SRes = 100; % search resolution
varphi_range = linspace(-pi/2,pi/2,SRes);
theta_range = linspace(-pi/2,pi/2,SRes);
a_varphi_range = zeros(M,SRes,SRes); % [M,Azimuth,Elevation]
% obtain the array response vectors for all azimuth-elevation pairs
parfor i = 1:SRes
    a_varphi_range(:,:,i) = ...
    UPA_Evaluate(lambda,M_V,M_H,varphi_range,repelem(theta_range(i),1,SRes),elementspacing,elementspacing);
end

idx = zeros(nbrOfAngleRealizations,nbrOfNoiseRealizations); 
for n1 = 1:nbrOfAngleRealizations
    disp(n1);

    %Select the parameters to be estimated 
    %channel g (UE to RIS)
    varphi_UE = azimuth(n1);
    theta_UE = elevation(n1);
    var_amp_g= 1;
    g = sqrt(var_amp_g) * Cph(n1) * UPA_Evaluate(lambda,M_V,M_H,varphi_UE,theta_UE,elementspacing,elementspacing);
    %channel d (UE to BS)   
    var_pow_d= 10;
    d = sqrt(var_pow_d/2) * (randn + 1i*randn); % CN(0,10)

    %Compute the exact capacity for the system Eq. (3)
    capacity(n1) = log2(1+SNR_data*(sum(abs(Dh*g)) + abs(d))^2);
    
    for n2 = 1:nbrOfNoiseRealizations
        if mod(n1,Prep) == 1
            %Select which two RIS configurations from the grid of beams to start with
            utilize = false(M,1);
            if n1 == 1
                % For the first time we initialize randomly 
                utilize(round(M/3)) = true;
                utilize(round(2*M/3)) = true;
                Plim = M; % number of pilots for intial estimation
            else
                % Start the estimation with using the previous best RIS
                % config to track the channel
                utilize(bestInit1) = true;
                idx(n1,n2) = randi(M-M_H)+M_H;
                while idx(n1,n2) == bestInit1
                    idx(n1,n2) = randi(M-M_H)+M_H;
                end
                utilize(idx(n1,n2)) = true;
                Plim = 6; % number of pilots for tracking
            end

            %Define the initial transmission setup
            RIS_directions = beamAngles(utilize,:);
            L = length(RIS_directions); 
            RISconfigs = Dh_angles*UPA_Evaluate(lambda,M_V,M_H,RIS_directions(:,2),RIS_directions(:,1),elementspacing,elementspacing);
            B = RISconfigs';

            %Generate the noise
            noise = (randn(M,1)+1i*randn(M,1))/sqrt(2);


            %Go through iterations by adding extra RIS configurations in the estimation
            for itr = 1:Plim-1
            %Generate the received signal
                y =  sqrt(SNR_pilot)*(B*Dh*g + d) + noise(1:itr+1,1);

                %Compute the ML utility function for all potential angle pairs
                utility_num = zeros(SRes,SRes); %numerator of the utility function
                utility_den = zeros(SRes,SRes); %denominator of the utility function
             
                % Each row is for one elevation angle
                for i = 1:SRes
                    utility_num(i,:) = abs(y' * (eye(itr+1) - (itr+1)^-1 * ones(itr+1,itr+1))* ...
                        B*Dh*a_varphi_range(:,:,i)).^2;
                    utility_den(i,:) = sum(abs(B*Dh*a_varphi_range(:,:,i)).^2,1) - (itr+1)^-1 * ...
                        abs(ones(1,itr+1)*B*Dh*a_varphi_range(:,:,i)).^2;
                end
           
             
                utilityfunction = utility_num ./ utility_den;

                %Extract the angle estimate
                [~,maxind] = max(utilityfunction,[],'all');
                [Elidx,Azidx] = ind2sub([SRes,SRes],maxind);

                %Estimate g
                var_amp_g_num = abs(y' * (eye(itr+1) - (itr+1)^-1 * ones(itr+1,itr+1))* B*Dh*...
                 UPA_Evaluate(lambda,M_V,M_H,varphi_range(Azidx),theta_range(Elidx),elementspacing,elementspacing))^2;
                var_amp_g_den = SNR_pilot * (sum(abs(B*Dh*UPA_Evaluate(lambda,M_V,M_H,varphi_range(Azidx),theta_range(Elidx),elementspacing,elementspacing)).^2,1) ...
                    - (itr+1)^-1 * abs(ones(1,itr+1)*B*Dh*UPA_Evaluate(lambda,M_V,M_H,varphi_range(Azidx),theta_range(Elidx),elementspacing,elementspacing))^2)^2;
                var_amp_g_est = var_amp_g_num/var_amp_g_den;
                var_phas_g_est = - angle (y' * (eye(itr+1) - (itr+1)^-1 * ones(itr+1,itr+1))* B*Dh*UPA_Evaluate(lambda,M_V,M_H,varphi_range(Azidx),theta_range(Elidx),elementspacing,elementspacing));
                g_est = sqrt(var_amp_g_est) * exp(1i*var_phas_g_est) * UPA_Evaluate(lambda,M_V,M_H,varphi_range(Azidx),theta_range(Elidx),elementspacing,elementspacing);

                %Estimate d
                var_amp_d_est = (SNR_pilot)^-1 * (itr+1)^-2 * abs (ones(1,itr+1) * (y - sqrt(SNR_pilot)*B*Dh*g_est))^2;
                var_phas_d_est = angle (ones(1,itr+1) * (y - sqrt(SNR_pilot)*B*Dh*g_est));
                d_est = sqrt(var_amp_d_est) * exp(1i*var_phas_d_est);
             
             
                %Estimate the RIS configuration that (approximately) maximizes the SNR
                RISconfig = angle(Dh*g_est)-var_phas_d_est;
     

                %Compute the corresponding achievable rate
                rate_proposed(itr,n1,n2) = log2(1+SNR_data*abs(exp(-1i*RISconfig).'*Dh*g + d)^2);


                %Find an extra RIS configuration to use for pilot transmission
                if itr < Plim-1

                    %Find which angles in the grid-of-beams haven't been used
                    unusedAngles = beamAngles(utilize==false,:);

                    %Guess what the channel would be with the different beams
                    guessOnAngles = Dh*UPA_Evaluate(lambda,M_V,M_H,unusedAngles(:,2),unusedAngles(:,1),elementspacing,elementspacing);
               


                    %Find which of the guessed channels matches best with the
                    %currently best RIS configuration

                    closestBeam = abs(exp(-1i*RISconfig).'*guessOnAngles); 
                    [~,bestUnusedBeam] = max(closestBeam);


                    %Add a pilot transmission using the new RIS configuration
                    [~,newAngle] = ismember(unusedAngles(bestUnusedBeam,:),beamAngles,"rows");
                    utilize(newAngle) = true;
                    RIS_directions = beamAngles(utilize,:);
                    RISconfigs = Dh_angles*UPA_Evaluate(lambda,M_V,M_H,RIS_directions(:,2),RIS_directions(:,1),...
                        elementspacing,elementspacing);
                    B = RISconfigs';
                else
                    closestBeam = abs(exp(-1i*RISconfig).'*Dh_angles*beamresponses); 
                    [~,bestInit1] = max(closestBeam);
                end
            end
        else
            rate_proposed(itr,n1,n2) = log2(1+SNR_data*abs(exp(-1i*RISconfig).'*Dh*g + d)^2);
        end
    end
end
Rred = mean(rate_proposed,3); % redundant variable to store rate values for plot
R = zeros(1,RWL);
R(1:Prep) = Rred(M-1,1:Prep);
R(Prep+1:end) = Rred(Plim-1,Prep+1:end); % collection of rate for all time
set(groot,'defaultAxesTickLabelInterpreter','latex');  
figure; 
plot(capacity,'-k','LineWidth',5); hold on; % plot upperbound
plot(R,'LineWidth',2); % plot the achievable rate
% plot the updating points
plot(1:Prep:nbrOfAngleRealizations,R(1:Prep:nbrOfAngleRealizations),'Marker','o','LineStyle','none','MarkerSize',15,'LineWidth',4);
xticks(0:100:1000); ylim([0,9.2]); grid on; box on;
xticklabels({'0','20','40','60','80','100','120','140','160','180','200'})
legend({'Perfect CSI','1\textsuperscript{st} Policy (1 second)','2\textsuperscript{nd} Policy (10 seconds)','Update points'},'Fontsize',20,'Interpreter','latex');
ylabel('Spectral Efficiency (bit/s/Hz)','Interpreter','latex');
xlabel('Time (s)','Interpreter','latex');