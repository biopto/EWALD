%% Tomographic Reconstruction with Gerchberg-Papoulis algorithm + Finite Object Support
% Authors: Piotr Makowski, Wojciech Krauze, Paweł Ossowski
% Warsaw University of Technology, 2024

clear; clc;
close all
[path, ~, ~] = fileparts(mfilename('fullpath'));
addpath(genpath(path));

%%%%%%%%% START OF PARAMETERS SECTION:

%% Reconstruction Presets
DI; ... % direct inversion (no iterations)
% GP; ... % iterations
% GPSC; ... % iterations + finite object support

%% Optimization parameters
% After filling all projections: there is no useful information in 'z' far from center, so we can crop it by downsampling Kspace before IFT - must be power of 2
ROI_crop_z = 4; % 1 / 2 / 4 (1 - no crop)
% Resolution in 'z' is 4x worse then 'xy' so before IFT, Kspace can be cropped in 'z' to 0.25 of its original size
limit_resolution_z = 0.25; % 0.25 / 0.5 / 1 (1 - no crop)
sinogram_subset_factor = 1;     % default = 1     % take every Nth projection from sinogram
projection_crop_factor = 1;     % default = 1;    % crop factor (<=1)
resample_projections   = false; % default = false % resample projections to minimal safe resolution

%% Save reconstruction?
save_reconstruction    = 1;     % 0 - don't save; 1 - save RI; 2 - save RI and absorption (complex RI)

%% Plots
% show selected plots:
% '0' : off
% '1' : show final reconstruction
% '2' : show spectrum of each processed projection
% '3' : show fully filled K-space, Direct Inversion result and object support (if applies)
% '4' : show reconstruction (XZ) progress in GP and 1D plot of relative progress of GP
% '5' : show selected plots of reconstruction output and save them to .png (Plots.m script)
% 'D' : debug mode - detailed plots not shown in other modes:
%                  - show the process of filling K-space with projections
% 'DD': extended debug mode:
%                  - show the process of filling K-space with projections
%                  - show the process of filling K-space with projections- also show spectrum of projection with dimmed out frequencies aroud NA
plots = '1';

%%%%%%%%% END OF PARAMETER SECTION

%% Load sinogram and system basic parameters
loadfiles;

%% Assume defaults for unspecified things
defaults;

%% Preparing data for reconstruction
[SINOamp_reduced,SINOph_reduced, sino_params, dx, ...
 projection_padding_xy, N_projection_padded, ...
 x_tv, n_obj_res, plots, Fpmask, projection_downsample_factor, N_SINO_y] = ...
	Preprocess(SINOamp, SINOph, geometry, sino_params, resample_projections, dx, projection_crop_factor, Fpmask,...
                n_obj, thetay, sinogram_subset_factor, n_imm, Masking, N_CP, nCPi,...
                projection_padding_xy, plots, do_NNC);

%% Calculating reconstruction   
[RECON, dx, dz, nGPi, N_Kspace_xy_padded, KO,KOi, RMAEtab,RRMSEtab,RMADtab,RRMSDtab] = ...
	FDT(SINOamp_reduced,SINOph_reduced, sino_params, thetay, ... % sinogram (including index of OCT proj)
		n_imm,dx, ... % optical system (ODT and OCT)
		geometry,Approx,interpFp,Ramp,do_NNC, ... % solver approximations
		projection_padding_xy, Kspace_padding, N_projection_padded, Kspace_oversampling_z, ROI_crop_z, limit_resolution_z, ... % Fourier space sampling
		nGPi,epsi,relaxGP,relaxM, ... % Gerchberg-Papoulis iterationsnGPi
		x_tv, n_obj_res, ... % Gerchberg-Papoulis object constraints/reference
		plots, Fpmask);

%% Display and save
time=toc;
sim_time = show_time(time);

if save_reconstruction
    save_recon;
end

if contains(plots,'1'); vis(real(RECON), [min(real(RECON(:))), max(real(RECON(:)))]); colormap jet; end
if contains(plots,'5'); Plots; end
