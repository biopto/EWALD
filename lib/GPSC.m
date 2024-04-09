% reconstruction parameters
% GPSC iterations

%% Solver parameters
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nGPi = 30;%30 % max number of GP iterations with Masking==1
relaxGP = 1.0;%1.0 % relaxation parameter for Fourier data replenishment (LaR: and rec padding suppression)
epsi = 10*0.8e-8;%10*0.8e-8(good saturation,<=nGPi), NaN(iters=nGPi) % epsilon value for auto-stop condition, terminate quicker with Masking

Masking = true;% toggle spatial support in GP iterations
relaxM = 1.0;%1.0 % relaxation parameter for spatial support
nCPi = 0;          % set 0 for auto % number of Chambolle-Pock iterations
N_CP = 250;         % size of Chambolle-Pock reconstruction from which object support is generated (250 is optimal)