close all; clc;
clear;
freq = 28e9; % Central frequency
lambda = physconst('LightSpeed') / freq; % Wavelength
SRes = 100; % search resolution
%UPA Element configuration
M_H = 16; M_V = 16; M = M_H*M_V;
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

nbrOfAngleRealizations = 100;
nbrOfNoiseRealizations = 10;


%Save the rates achieved at different iterations of the algorithm
capacity = zeros(1,nbrOfAngleRealizations);
rate_proposed = zeros(M-1,nbrOfAngleRealizations,nbrOfNoiseRealizations);
rate_LS = zeros(M-1,nbrOfAngleRealizations,nbrOfNoiseRealizations);


%Create a uniform grid of beams (like a DFT matrix) to be used at RIS
ElAngles = asin((-M_V/2:1:M_V/2-1)*2/M_V);
AzAngles = asin((-M_H/2:1:M_H/2-1)*2/M_H);
beamAngles = zeros(M,2);

for i = 1:length(ElAngles)
    figure(1);
    plot(rad2deg(AzAngles),rad2deg(repelem(ElAngles(i),1,length(AzAngles))),'*');
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
    a_varphi_range = zeros(M,SRes,SRes);
    for i = 1:SRes
        a_varphi_range(:,:,i) = ...
            UPA_Evaluate(lambda,M_V,M_H,varphi_range,repelem(theta_range(i),1,SRes),elementspacing,elementspacing);
    end

    %Compute the exact capacity for the system Eq. (3)
    capacity(n1) = log2(1+SNR_data*sum(abs(Dh*g))^2);

    for n2 = 1:nbrOfNoiseRealizations
        L = 2;
        %Select which two RIS configurations from the grid of beams to start with
        utilize = false(1,M);
        RISconfigs = zeros(M,L);
        for l = 1:L
            redunt = round(l*M/3);
            utilize(redunt) = true;
            AzSel = AzAngles(mod(redunt-1,M_H)+1); ElSel = ElAngles(floor(redunt/M_H)+1);
            figure(1);
            plot(rad2deg(AzSel),rad2deg(ElSel),'o','MarkerSize',15);
            RISconfigs(:,l) = UPA_Evaluate(lambda,M_V,M_H,AzSel,ElSel,elementspacing,elementspacing);
        end

        %Define the initial transmission setup
        RISconfigs = Dh_angles*RISconfigs; % ????????
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
            utilityfunction = abs(y'*B*Dh*a_varphi_range).^2./sum(abs(B*Dh*a_varphi_range).^2,1);


            %Extract the angle estimate
            [~,maxind] = max(utilityfunction);


            %Estimate the RIS configuration that (approximately) maximizes the SNR
            RISconfig = angle(Dh*arrayresponse(varphi_range(maxind),M));

            %Compute the corresponding achievable rate
            rate_proposed(itr,n1,n2) = log2(1+SNR_data*abs(exp(-1i*RISconfig).'*Dh*g).^2);


            %Find an extra RIS configuration to use for pilot transmission
            if itr < M-1

                %Find which angles in the grid-of-beams that haven't been used
                unusedAngles = beamAngles(utilize==false);

                %Guess what the channel would be with the different beams
                guessOnAngles = Dh*arrayresponse(unusedAngles,M);


                %Find which of the guessed channels that matches best with the
                %currently best RIS configuration

                closestBeam = abs(exp(-1i*RISconfig).'*guessOnAngles); %Emil's suggestion
                [~,bestUnusedBeam] = max(closestBeam);


                %Add a pilot transmission using the new RIS configuration
                newAngle = find(unusedAngles(bestUnusedBeam) == beamAngles);
                utilize(newAngle) = true;
                RIS_directions = beamAngles(utilize);
                RISconfigs = Dh_angles*arrayresponse(RIS_directions,M);
                B = RISconfigs';
            end

        end

        %LS estimation
        DFT = fft(eye(M));
        randomOrdering = randperm(M);

        %Go through iterations by adding extra RIS configurations in the estimation
        for itr = 1:M-1

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