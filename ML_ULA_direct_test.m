%Testing the case with direct path 
%close all;
clear;

%Number of horizontal elements in a ULA
M = 40;
elementspacing = 1/2; %In wavelengths

%Set the SNR
SNR_dB = 0;
SNR = db2pow(SNR_dB);

%Select angle to the base station (known value)
varphi_BS = -pi/6;

%Select angle to the user (to be estimated)
varphi_UE = pi/3;

%Select the amplitude and phase of g and d
var_amp_g=1;
var_amp_d=0.1;
var_phas_g=pi/3;
var_phas_d=pi/3;

%Define the array response vector
arrayresponse = @(phi,M) exp(-1i*2*pi*elementspacing*(0:M-1)'*sin(phi));

%Generate channels
h = arrayresponse(varphi_BS,M); 
Dh = diag(h); Dh_angles = diag(h./abs(h));
g = sqrt(var_amp_g) * exp (1i*var_phas_g) * arrayresponse(varphi_UE,M);
d = sqrt(var_amp_d) * exp (1i*var_phas_d);


%Create a uniform grid of beams (like a DFT matrix) to be used at RIS
beamAngles = asin((-M/2:1:M/2)*2/M);

%Select which two RIS configurations from the grid of beams to start with
utilize = false(1,M);
utilize(round(M/3)) = true;
utilize(round(2*M/3)) = true;


%Define the initial transmission setup
RIS_directions = beamAngles(utilize);
RISconfigs = Dh_angles*arrayresponse(RIS_directions,M);
B = RISconfigs';


%% Define a fine grid of angle directions to analyze when searching for angle of arrival
varphi_range = linspace(-pi/2,pi/2,1000);

%Compute the exact capacity for the system
capacity = log2(1+SNR*(sum(abs(Dh*g)) + abs(d))^2);

%Save the rates achieved at different iterations of the algorithm
rate = zeros(M-1,1);

%Generate the noise
noise = (randn(M,1)+1i*randn(M,1))/sqrt(2);


%Go through iterations by adding extra RIS configurations in the estimation
for itr = 1:M-1

    %Generate the received signal
    y = sqrt(SNR)*(B*Dh*g + d*ones(itr+1,1)) + noise(1:itr+1,1);

    %Compute the ML utility function for all potential angles
    utility_num = zeros(length(varphi_range),1);
    utility_den = zeros(length(varphi_range),1);

    for i = 1:length(varphi_range)
        utility_num(i) = abs(y' * (eye(itr+1) - (itr+1)^-1 * ones(itr+1,itr+1))* B*Dh*arrayresponse(varphi_range(i),M))^2;
        utility_den(i) = sum(abs(B*Dh*arrayresponse(varphi_range(i),M)).^2,1) - (itr+1)^-1 * abs(ones(1,itr+1)*B*Dh*arrayresponse(varphi_range(i),M))^2;
  
    end
      utilityfunction = utility_num ./ utility_den;

    %Extract the angle estimate
    [~,maxind] = max(utilityfunction);
    arrayresponse(varphi_range(maxind),M);
    
    %Estimate the amplitude and phase of g
    var_amp_g_num = abs(y' * (eye(itr+1) - (itr+1)^-1 * ones(itr+1,itr+1))* B*Dh*arrayresponse(varphi_range(maxind),M))^2;
    var_amp_g_den = SNR * (norm(B*Dh*arrayresponse(varphi_range(maxind),M))^2 - (itr+1)^-1 * abs(ones(1,itr+1)*B*Dh*arrayresponse(varphi_range(maxind),M))^2)^2;
    var_amp_g_est = var_amp_g_num/var_amp_g_den;
    var_phas_g_est = - angle (y' * (eye(itr+1) - (itr+1)^-1 * ones(itr+1,itr+1))* B*Dh*arrayresponse(varphi_range(maxind),M));
    g_est = sqrt(var_amp_g_est) * exp(1i*var_phas_g_est) * arrayresponse(varphi_range(maxind),M);
    
    %Estimate the amplitude and phase of d
    z = (y - sqrt(SNR)*B*Dh*g_est);
    var_amp_d_est = (SNR)^-1 * (itr+1)^-2 * abs (ones(1,itr+1) * (y - sqrt(SNR)*B*Dh*g_est))^2;
    var_phas_d_est = angle (ones(1,itr+1) * (y - sqrt(SNR)*B*Dh*g_est));
    d_est = sqrt(var_amp_d_est) * exp(1i*var_phas_d_est);

    %Estimate the RIS configuration that (approximately) maximizes the SNR
    RISconfig = angle(Dh*g_est) - var_phas_d_est * ones (M,1) ;

    %Compute the corresponding achievable rate
    rate(itr) = log2(1+SNR*abs(exp(-1i*RISconfig).'*Dh*g + d)^2);

   %Find an extra RIS configuration to use for pilot transmission
   if itr < M-1

   %Find the angles in the grid-of-beams that haven't been used
   unusedAngles = beamAngles(utilize==false);

   %Guess what the channel would be with the different beams
   guessOnAngles = Dh*exp(1i * var_phas_g_est)*arrayresponse(unusedAngles,M);

   %Find which of the guessed channels matches best with the currently best RIS configuration
   closestBeam = abs(exp(-1i*RISconfig).'*guessOnAngles + d_est);
   [~,bestUnusedBeam] = max(closestBeam);

   %Add a pilot transmission using the new RIS configuration
   newAngle = find(unusedAngles(bestUnusedBeam) == beamAngles);
   utilize(newAngle) = true;
   RIS_directions = beamAngles(utilize);
   phase_shifts = angle (Dh_angles*exp(1i * var_phas_g_est)*arrayresponse(RIS_directions,M)) - var_phas_d_est;
   RISconfigs = exp (1i * phase_shifts);
   B = RISconfigs';
   end

end


figure;
hold on; box on; grid on;
plot(1:M-1,capacity*ones(M-1,1),'r:','LineWidth',2)
plot(1:M-1,rate,'k-.','LineWidth',2)
ax = gca;
xlabel('Iterations','Interpreter','latex');
ylabel('Rate','Interpreter','latex');
set(gca,'fontsize',16);






