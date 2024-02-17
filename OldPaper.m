close all;
clear;

%Number of horizontal elements in a ULA
M = 64;
d_H = 1/2; %In wavelengths
Plim = 16;
LS = false;
%Set the SNR
SNRdB_pilot = 10;
SNR_pilot = db2pow(SNRdB_pilot);

SNRdB_data = 0;
SNR_data = db2pow(SNRdB_data);

%Select angle to the base station (known value)
varphi_BS = -pi/6;

%Define the array response vector
arrayresponse = @(phi,M) exp(-1i*2*pi*d_H*(0:M-1)'*sin(phi));

%Generate channels
h = arrayresponse(varphi_BS,M);
Dh = diag(h);
Dh_angles = diag(h./abs(h));

nbrOfAngleRealizations = 10;
nbrOfNoiseRealizations = 10;


%Save the rates achieved at different iterations of the algorithm
capacity = zeros(1,nbrOfAngleRealizations);
rate_proposed = zeros(Plim-1,nbrOfAngleRealizations,nbrOfNoiseRealizations);
rate_LS = zeros(Plim-1,nbrOfAngleRealizations,nbrOfNoiseRealizations);


%Create a uniform grid of beams (like a DFT matrix) to be used at RIS
beamAngles = asin((-M/2:1:M/2)*2/M);

% Define a fine grid of angle directions to analyze when searching for angle of arrival
varphi_range = linspace(-pi/2,pi/2,1000);
a_varphi_range = arrayresponse(varphi_range,M); 
utilityfunction = zeros(1000,Plim);    
for n1 = 1:nbrOfAngleRealizations
    n1

    %Select angle to the user (to be estimated)
    varphi_UE = rand(1)*2*pi/3-pi/3;
    g = arrayresponse(varphi_UE,M);


    %Compute the exact capacity for the system
    capacity(n1) = log2(1+SNR_data*sum(abs(Dh*g))^2);

    for n2 = 1:nbrOfNoiseRealizations

        %Select which two RIS configurations from the grid of beams to start with
        utilize = false(1,M);
        utilize(round(M/3)) = true;
        utilize(round(2*M/3)) = true;


        %Define the initial transmission setup
        RIS_directions = beamAngles(utilize);
        L = length(RIS_directions);
        RISconfigs = Dh_angles*arrayresponse(RIS_directions,M);
        B = RISconfigs';

        %Generate the noise
        noise = (randn(M,1)+1i*randn(M,1))/sqrt(2);


        %Go through iterations by adding extra RIS configurations in the estimation
        for itr = 1:Plim-1

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
            utilityfunction(:,itr) = transpose(abs(y'*B*Dh*a_varphi_range).^2./sum(abs(B*Dh*a_varphi_range).^2,1));

            %Extract the angle estimate
            [~,maxind] = max(utilityfunction(:,itr));


            %Estimate the RIS configuration that (approximately) maximizes the SNR
            RISconfig = angle(Dh*arrayresponse(varphi_range(maxind),M));

            %Compute the corresponding achievable rate
            rate_proposed(itr,n1,n2) = log2(1+SNR_data*abs(exp(-1i*RISconfig).'*Dh*g).^2);


            %Find an extra RIS configuration to use for pilot transmission
            if itr < Plim-1

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

        if LS
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

end


set(groot,'defaultAxesTickLabelInterpreter','latex');  

figure;
hold on; box on; grid on;
plot(2:Plim,mean(capacity)*ones(Plim-1,1),'r:','LineWidth',2)
plot(2:Plim,mean(mean(rate_proposed,3),2),'k-','LineWidth',2)
plot(2:Plim,mean(mean(rate_LS,3),2),'b-.','LineWidth',2)
ax = gca;
xlabel('Number of pilot transmissions','Interpreter','latex');
ylabel('Average rate (bits/s/Hz)','Interpreter','latex');
legend({'Perfect CSI','Proposed ML estimator','Least-squares estimator'},'Interpreter','latex','Location','SouthEast');
set(gca,'fontsize',16);
ylim([0 ceil(max(capacity))+0.5]);





