close all; clc; clear; 
freq = 28e9; % Central frequency
lambda = physconst('LightSpeed') / freq; % Wavelength

% Mobility Config 
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

SRes = 100; % search resolution
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

nbrOfAngleRealizations = 1000;
nbrOfNoiseRealizations = 1;


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

%plot the Configured angles grid
% for i = 1:length(beamAngles)
%     figure(1);
%     grid on;
%     plot(rad2deg(beamAngles(i,2)),rad2deg(beamAngles(i,1)),'*','MarkerSize',10,'Color','b');
%     hold on;
% end

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

    % Define a fine grid of angle directions to analyze when searching for angle of arrival
    varphi_range = linspace(-pi/2,pi/2,SRes);
    theta_range = linspace(-pi/2,pi/2,SRes);
    a_varphi_range = zeros(M,SRes,SRes); % [M,Azimuth,Elevation]

    % obtain the array response vectors for all azimuth-elevation pairs
    parfor i = 1:SRes
        a_varphi_range(:,:,i) = ...
            UPA_Evaluate(lambda,M_V,M_H,varphi_range,repelem(theta_range(i),1,SRes),elementspacing,elementspacing);
    end

    %Compute the exact capacity for the system Eq. (3)
    capacity(n1) = log2(1+SNR_data*(sum(abs(Dh*g)) + abs(d))^2);
    if mod(n1,Prep) == 1
        for n2 = 1:nbrOfNoiseRealizations

            %Select which two RIS configurations from the grid of beams to start with
            utilize = false(M,1);
            utilize(round(M/3)) = true;
            utilize(round(2*M/3)) = true;

            %Define the initial transmission setup
            RIS_directions = beamAngles(utilize,:);
            L = length(RIS_directions); 
            RISconfigs = Dh_angles*UPA_Evaluate(lambda,M_V,M_H,RIS_directions(:,2),RIS_directions(:,1),elementspacing,elementspacing);
            B = RISconfigs';

            %Generate the noise
            noise = (randn(M,1)+1i*randn(M,1))/sqrt(2);


            %Go through iterations by adding extra RIS configurations in the estimation
            for itr = 1:M-1
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
                if itr < M-1

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
                end
            end
        end

        %LS estimation
%         DFT = fft(eye(M+1));
%         randomOrdering = randperm(M+1);
        
%         %Go through iterations by adding extra RIS configurations in the estimation
%         for itr = 1:M-1
% 
%             B_LS = transpose(DFT(:,randomOrdering(1:itr+1)));
% 
%             v = [d;Dh*g]; % cascaded channel and direct path 
%             y = sqrt(SNR_pilot)*B_LS*v + noise(1:itr+1,1);
%              
%             %Compute LS estimate without parametrization
%             v_LS = pinv(B_LS)*y/sqrt(SNR_pilot);
% 
%             %Estimate the RIS configuration that (approximately) maximizes the SNR
%             RISconfig =  angle(v_LS(2:end)) - angle(v_LS(1));
% 
%             %Compute the corresponding achievable rate
%             rate_LS(itr,n1,n2) = log2(1+SNR_data*abs(exp(-1i*RISconfig).'*Dh*g+d).^2);
% 
%         end

    else
        rate_proposed(itr,n1,n2) = log2(1+SNR_data*abs(exp(-1i*RISconfig).'*Dh*g + d)^2);
        capacity(n1) = log2(1+SNR_data*(sum(abs(Dh*g)) + abs(d))^2);
    end
end

set(groot,'defaultAxesTickLabelInterpreter','latex');  
figure; 
plot(capacity,'--','LineWidth',2); hold on;
plot(rate_proposed(63,:),'LineWidth',2);
%xline(1:Prep:nbrOfAngleRealizations);
plot(1:Prep:nbrOfAngleRealizations,rate_proposed(63,1:Prep:nbrOfAngleRealizations),'Marker','o','LineStyle','none','MarkerSize',10,'LineWidth',2);
xticks(0:100:1000);
xticklabels({'0','20','40','60','80','100','120','140','160','180','200'})
legend({'Perfect CSI','1 seconds','10 seconds','Update '},'Fontsize',20,'Interpreter','latex');

figure;
hold on; box on; grid on;
plot(2:M,mean(capacity)*ones(M-1,1),'r:','LineWidth',2)
plot(2:M,mean(mean(rate_proposed,3),2),'k-','LineWidth',2)
plot(2:M,mean(mean(rate_LS,3),2),'b-.','LineWidth',2)
ax = gca;
xlabel('Number of pilot transmissions','Interpreter','latex');
ylabel('Average rate (bits/s/Hz)','Interpreter','latex');
legend({'Perfect CSI','Proposed ML estimator','Least-squares estimator'},'Interpreter','latex','Location','SouthEast');
set(gca,'fontsize',16);
ylim([0 ceil(max(capacity))+0.5]);
