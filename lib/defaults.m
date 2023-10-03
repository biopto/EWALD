if exist('n_obj', 'var')
	n_obj(isnan(n_obj)) = n_imm; % immerse the synthetic object
else
	n_obj = [];
end

if ~exist('do_NNC', 'var')
	warning('`do_NNC` missing in parameters. Assuming do_NNC=1')
	do_NNC = 1; % non-negativity constraint by default
end

if ~exist('geometry', 'var')
	warning('`geometry` missing in parameters. Assuming `fixed`')
	geometry = 'fixed'; % LAT configuration
end

if ~exist('dx', 'var') % dx and everything depending on M should be expressed in sample plane already
    dx = cam_pix*downsampling/M;
end

if ~exist('thetay', 'var') % angles for full-angle object-rotation mode around Y
	thetay = [];
else
	kx_red = []; ky_red = [];
end

if ~exist('Fpmask','var')
    Fpmask = [];
elseif ismatrix(Fpmask)
    Fpmask = repmat(Fpmask, [1 1 size(SINOph,3)]);
end

if Masking==0
    N_CP = [];
    nCPi = [];
end

% backward compatibility
% sino_params:
% 1-transmission, 2-reflection
if ~exist('sino_params','var')
    sino_params(1:2,:) = [rayXY]; % illumination directions of each projection
    sino_params(3,  :) = ones(1,size(rayXY,2)).*lambda; % wavelength of each projection
    sino_params(4,  :) = ones(1,size(rayXY,2)).*NA; % numerical aperture of each projection
    sino_params(5,  :) = ones(1,size(rayXY,2)); % default is transmission system for each projection
end

if ~exist('Kspace_padding','var')
    Kspace_padding = 'optimal'; % 'optimal', 'maxaperture' 
else
    if ~(strcmp(Kspace_padding,'optimal') || strcmp(Kspace_padding,'maxaperture'))
        error('Kspace_padding must be either "optimal" or "maxaperture"')
    end
end

tic;