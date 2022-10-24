close all; clc;
clear;
freq = 28e9; % Central frequency
lambda = physconst('LightSpeed') / freq; % Wavelength
SRes = 100; % search resolution
%UPA Element configuration
M_H = 8; M_V = 8; M = M_H*M_V;
elementspacing = 1/4; %In wavelengths

%Set the SNR
SNRdB_pilot = 10;
SNR_pilot = db2pow(SNRdB_pilot);

SNRdB_data = 0;
SNR_data = db2pow(SNRdB_data);

%Select angle to the base station (known value)
varphi_BS = -pi/6;
theta_BS = 0;

%Generate channels
h = UPA_Evaluate(lambda,M_V,M_H,varphi_BS,theta_BS,elementspacing,elementspacing);
Dh = diag(h);
Dh_angles = diag(h./abs(h));

nbrOfAngleRealizations = 2;
nbrOfNoiseRealizations = 10;


%Save the rates achieved at different iterations of the algorithm
capacity = zeros(1,nbrOfAngleRealizations);
rate_proposed = zeros(M-1,nbrOfAngleRealizations,nbrOfNoiseRealizations);
rate_LS = zeros(M-1,nbrOfAngleRealizations,nbrOfNoiseRealizations);


%Create a uniform grid of beams (like a DFT matrix) to be used at RIS
ElAngles = asin((-M_V/2:1:M_V/2-1)*2/M_V);
AzAngles = asin((-M_H/2:1:M_H/2-1)*2/M_H);
beamAngles = zeros(M,2); % Elevation-Azimuth pair
for i = 1:length(ElAngles)
    for j = 1:length(AzAngles)
        beamAngles(length(AzAngles)*(i-1)+j,:) = [ElAngles(i) AzAngles(j)];
    end
end

% plot the Configured angles grid
for i = 1:length(beamAngles)
    figure(1);
    grid on;
    plot(rad2deg(beamAngles(i,2)),rad2deg(beamAngles(i,1)),'*','MarkerSize',10,'Color','b');
    hold on;
end

for n1 = 1:nbrOfAngleRealizations
    disp(n1);
    %Select angle to the user (to be estimated)
    varphi_UE = rand(1)*2*pi/3-pi/3;
    theta_UE = -rand(1)*pi/2;
    g = UPA_Evaluate(lambda,M_V,M_H,varphi_UE,theta_UE,elementspacing,elementspacing);

    % Define a fine grid of angle directions to analyze when searching for angle of arrival
    varphi_range = linspace(-pi/2,pi/2,SRes);
    theta_range = linspace(-pi/2,0,SRes);
    a_varphi_range = zeros(M,SRes,SRes); % [M,Azimuth,Elevation]
    parfor i = 1:SRes
        a_varphi_range(:,:,i) = ...
            UPA_Evaluate(lambda,M_V,M_H,varphi_range,repelem(theta_range(i),1,SRes),elementspacing,elementspacing);
    end

    %Compute the exact capacity for the system Eq. (3)
    capacity(n1) = log2(1+SNR_data*sum(abs(Dh*g))^2);

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
            y = sqrt(SNR_pilot)*B*Dh*g + noise(1:itr+1,1);

            %Compute the ML utility function for all potential angles
%            utilityfunction = zeros(length(varphi_range),1);

%             for i = 1:length(varphi_range)
% 
%                 utilityfunction(i) = abs(y'*B*Dh*a_varphi_range(:,i))^2/norm(B*Dh*a_varphi_range(:,i))^2;
% 
%             end

            %Compute the ML utility function for all potential angles - This new version is 50 times faster!
            % Each row is for one elevation angle
            utilityfunction = zeros(SRes,SRes);
            parfor i = 1:SRes
                utilityfunction(i,:) = abs(y'*B*Dh*a_varphi_range(:,:,i)).^2./sum(abs(B*Dh*a_varphi_range(:,:,i)).^2,1);
            end
            %Extract the angle estimate
            [~,maxind] = max(utilityfunction,[],'all');
            [Elidx,Azidx] = ind2sub([SRes,SRes],maxind);

            %Estimate the RIS configuration that (approximately) maximizes the SNR
            RISconfig = angle(Dh*UPA_Evaluate(lambda,M_V,M_H,varphi_range(Azidx),theta_range(Elidx),elementspacing,elementspacing));

            %Compute the corresponding achievable rate
            rate_proposed(itr,n1,n2) = log2(1+SNR_data*abs(exp(-1i*RISconfig).'*Dh*g).^2);


            %Find an extra RIS configuration to use for pilot transmission
            if itr < M-1

                %Find which angles in the grid-of-beams that haven't been used
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
                RISconfigs = Dh_angles*UPA_Evaluate(lambda,M_V,M_H,RIS_directions(:,2),RIS_directions(:,1),elementspacing,elementspacing);
                B = RISconfigs';
            end

        end

        %LS estimation
        DFT = fft(eye(M));
        randomOrdering = randperm(M);

        %Go through iterations by adding extra RIS configurations in the estimation
        parfor itr = 1:M-1

            B_LS = transpose(DFT(:,randomOrdering(1:itr+1)));

            %Generate the received signal
            y = sqrt(SNR_pilot)*B_LS*Dh*g + noise(1:itr+1,1);

            %Compute LS estimate without parametrization
            g_LS = (Dh\(pinv(B_LS)*y))/sqrt(SNR_pilot);

            %Estimate the RIS configuration that (approximately) maximizes the SNR
            RISconfig = angle(Dh*g_LS);

            %Compute the corresponding achievable rate
            rate_LS(itr,n1,n2) = log2(1+SNR_data*abs(exp(-1i*RISconfig).'*Dh*g).^2);

        end

    end

end


set(groot,'defaultAxesTickLabelInterpreter','latex');  

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