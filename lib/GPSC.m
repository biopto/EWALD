% reconstruction parameters
% GPSC iterations

%% Solver parameters
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Approx = 'Rytov';% Born,Rytov,PhaseRay(non-diffractive) % weak scattering approximation
Ramp = false;%0 % Ram-Lak filter on/off (FAT, nGPi==0)

nGPi = 30;%30 % max number of GP iterations with Masking==1
relaxGP = 1.0;%1.0 % relaxation parameter for Fourier data replenishment (LaR: and rec padding suppression)
epsi = 10*0.8e-8;%10*0.8e-8(good saturation,<=nGPi), NaN(iters=nGPi) % epsilon value for auto-stop condition, terminate quicker with Masking

% First, upsample Kspace in 'z' direction to guarantee precise positioning of projections on Ewald spheres
Kspace_oversampling_z = 1.000;%4.0,>=2.0(MAP),1.375(GRID-1D) % 1-D Fourier oversampling along Z axis (efficient with MAP)
% After filling all projections: there is no useful information in 'z' far from center, so we can crop it by downsampling Kspace before IFT - must be power of 2
ROI_crop_z = 4; % 1 / 2 / 4 (1 - no crop)
% Finally, resolution in 'z' is 4x worse then 'xy' so before IFT, Kspace can be cropped in 'z' to 0.25 of its original size
limit_resolution_z = 0.25; % 0.25 / 0.5 / 1 (1 - no drop)

projection_padding_xy = 1.000;%1.0(alfaz>1.0),1.375(alfaz==1:GRID-3D) % projection padding factor for 3-D Fourier oversampling
% Interpolation method for Fp(NA) resampling before mapping onto Ewald mesh
% NOTE: with N_Kspace_z_padded_upsampled > N_Kspace_xy_padded only 'none' is useful!
interpFp = 'none'; % none(no resampling), linc(phase tilt), {nearest,linear,cubic,spline, ...}(interp2 methods)

Masking = true;% toggle spatial support in GP iterations
relaxM = 1.0;%1.0 % relaxation parameter for spatial support
nCPi = 0;          % set 0 for auto % number of Chambolle-Pock iterations
N_CP = 250;         % size of Chambolle-Pock reconstruction from which object support is generated (250 is optimal)