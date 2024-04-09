% reconstruction parameters
% GP iterations

%% Solver parameters
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nGPi = 50;%200 % max number of Gerchberg-Papoulis iterations with Masking==0
relaxGP = 1.0;%1.0 % relaxation parameter for Fourier data replenishment (LaR: and rec padding suppression)
epsi = 0.8e-8;%0.8e-8(good saturation,<=nGPi), NaN(iters=nGPi) % epsilon value for auto-stop condition

Masking = false;% toggle spatial support in GP iterations
relaxM = [];%1.0 % relaxation parameter for spatial support