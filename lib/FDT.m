% -------------------------------------------------------------------------
% Reconstruction engine for 3D Diffraction Tomography in transmission and reflection mode
% Copyright: 
% 2015-2018 Piotr L. Makowski, 
% 2019-2021 Wojciech Krauze, Paweł Ossowski
% -------------------------------------------------------------------------
function [RECON,dxo_xy, nGPi, N_Kspace_xy, KO,KOi, RMAEtab,RRMSEtab,RMADtab,RRMSDtab] = ...
			FDT(SINOamp,SINOph, sino_params, thetay, ... % sinogram
						n_imm,dx_projection, ... % optical system
						geometry,Approx,interpFp,Ramp,do_NNC, ... % solver approximations
						projection_padding_xy, Kspace_padding, N_projection, Kspace_oversampling_z, ROI_crop_z, limit_resolution_z,... % Fourier space sampling 
						nGPi,epsi,relaxGP,relaxM, ... % Gerchberg-Papoulis iterations
						x_tv, obj_type, n_obj, ... % Gerchberg-Papoulis object constraints/reference
						plots, Fpmask)

% (even sizes everywhere)
%%%%%%%%%%%% INPUT
% SINOamp                      - (x,y,proj) stack of absolute value of complex projections  
% SINOph                       - (x,y,proj) stack of unwrapped phase of complex projections
% sino_params                  - parameters of each projection in the sinogram
% sino_params(1,:)             - x-coordinates of illumination vectors
% sino_params(2,:)             - y-coordinated of illumination vectors    
%                              - limited-angle illumination scenario given by components of illuminating
%                                "wave vectors" |k|=n_imm/lambda, coresponding to projections in SINO(x,y,:);
%                                ignored if ~isempty(thetay);
% sino_params(3,:)             - wavelengths of illumination beam in vacuum
% sino_params(4,:)             - numerical apertures within object immersion liquid; NA = n_imm*sin(alpha_max)
% sino_params(5,:)             - measurement mode: 1-transmission, 2-reflection
% thetay                       - full-angle illumination scenario given by array of projection angles measured
%                                around Y axis, coresponding to projections in SINO(x,y,:);
%                                ~isempty(thetay) overrides sino_params(1,:) and sino_params(2,:)
% n_imm                        - refractive index of object immersion medium
% dx_projection                - projection sample size (CCD_pixel*sino_downsampling/system_magnification)
% projection_padding_xy		   - K space oversampling factor defining reconstruction space padding
% Kspace_padding               - parameter describing how the Kspace will be padded:
%                                "optimal" - the Kspace is padded optimally to hold all information 
%                                "maxaperture" - the Kspace is set to maximum size for the given NA and lambda
% N_projection_padded          - {>=Nx} projection size after padding, improves accuracy for Interp='MAP' (inefficient)
% Kspace_oversampling_z        - parameter describing oversampling of Kspace in 'z' direction for better mapping of data onto Ewald sphere
% geometry                     - capture system configuration (detector orientation):
%                                'facing' = rotating object  (kxp,kyp or thetay: detector plane follows projection wavefront)
%                                'fixed' = rotating illumination (kxp,kyp: detector plane fixed in XY)
% Approx                       - {'Born'|'Rytov'} weak scattering approximation
% Ramp                         - ~isempty(thetay): use Ram-Lak filter on every projection (along X dimension)
%				                 instead of averaging repeated voxels to mimic FBP or FBPR;
%				                 isempty(thetay): use 2D ramp filter optimized for conical illumination scenario
%                                NOTE: Ramp is useful only with direct inversion (nGPi==0)
% nGPi                         - number of data replenishment iterations (GPi) on K space
% do_NNC                       - apply non-negativity constraint (NNC):
%						         +1: normal NNC (n_obj>n_imm)
%						          0: disabled
%						       -1: inversed NNC (n_obj<n_imm)
% epsi                         - {<1|NaN} epsilon value for GP auto-stop condition
% init_obj                     - {(N_Kspace_xy_padded, N_Kspace_xy_padded, N_Kspace_z_padded_upsampled)|(Nx,Nx,Nx)|empty:[]} 
%                                3D matrix (x,y,z) to initialize n_rec before GP iterations
% x_tv                         - {(Nx,Nx,Nx)|empty:[]} 3D matrix (x,y,z) of the object spatial support, if known
% n_obj                        - {(Nx,Nx,Nx)|empty:[]} 3D matrix (x,y,z) of the object refractive index, if known
% relaxGP                      - (0..1) relaxation coefficient for Fourier data replenishment:
%						         beta[i+1]=relaxGP*beta[i], where beta[1]=1 is the starting replenishment factor
% relaxM                       - (0..1) relaxation coefficient for masking weight:
%						         beta[i+1]=relaxM*beta[i], where beta[1]=1 is the starting replenishment factor
% plots                        - toggle plots
% relaxLaR                     - relaxed suppression of reconstruction padding  from Fourier oversampling (<=1, 0:disable)
% Fpmask                       - stack of masks (or a single mask) for spectral filtering of projections under Rytov approximation;
%				                 use this  instead of passing down-pass-filtered projections to avoid wrapping of Rytov spectra tails;               
%				                 also enables GP to reconstruct corresponding unknown regions in 3D spectrum
%
%%%%%%%%%%%% OUTPUT
% n_rec                         - (x,y,z) reconstructed Re{complex_refractive_index} 
% dxo                           - object space sample size in XY
% nGPi                          - data replenishment iterations completed
% N_Kspace_xy_padded            - size of KO and n_rec arrays in XY
% KO                            - projection data in 3D Fourier domain
% KOi                           - 3D Fourier data after GP iterations
% RMAEtab                       - relative mean absolute error curve (NaN if isempty(n_obj) or nGPi==0)
% RRMSEtab                      - relative root mean square error curve (NaN if isempty(n_obj) or nGPi==0)
% RMADtab                       - relative mean absolute difference curve (NaN if nGPi==0)
% RRMSDtab                      - relative root mean square difference curve (NaN if nGPi==0)

%   * modified Gerchberg-Papoulis (GP) iterations:
%		- transparency constraint (TC): imag(n_rec)=0
%		- non-negativity constraint (NNC): real(n_rec)>=n_imm
%		- spatial support: n_rec(outside_object)=n_imm

if isempty(Fpmask); Fpmask0=[]; Fpmask=true(size(SINOamp(:,:,1))); else; Fpmask0=Fpmask; end

SINOamp = single(SINOamp);
SINOph = single(SINOph);
[Nx, Ny, nproj] = size(SINOph);
if Nx~=Ny
	error('Square projections required.');
elseif mod(Nx,2)
	SINOamp = SINOamp(1:end-1, 1:end-1, :); 
	SINOph = SINOph(1:end-1, 1:end-1, :); 
	[Nx, Ny, nproj] = size(SINOph);
end
% ---------------------------------------------------------------------------------------------------------
% apply ifftshift to projections before FFT to get smoother phase for Fourier mapping
objshift = true;% true(correct),false(projection spectra oscillating rapidly: mapping inaccurate)

qFp = 1;%1(none,spline),>1(others) % upsampling factor for Fp prior interpolation (integer)
% progress params
Nc = 120;% cropped cube size for incremental difference evaluation
% apply transparency constraint: imag(n_rec)=0
do_TC = true;%true
tukr = 0.25;%0.25  % Tukey window radius for projections (Approx=Born, interpFp=linc)
% relaxed suppression of reconstruction padding  from Fourier oversampling (<=1, 0:disable)
relaxLaR = 0;%1.00, 0.99(best with N_Kspace_z_padded_upsampled > N_Kspace_xy_padded), 0:disable !!!!
% ignore immersion volume in evaluation of quality measures (if n_obj given)
do_qmask = false; %false
% apply constraints to reconstruction after the final inverse FFT
do_constraint_rec = false;%true
% projection spectrum upsampling for interpolation (integer preferably)
if ~strcmp(interpFp,'none'); q = qFp;
else; q = 1; end
% data replenishment params
betaGP = 1.0;%1.0 % starting data replenishment weight (relaxed by relaxGP)
betaM = 1.0;%1.0 % starting mask weight (relaxed by relaxM)

% with newer versions of MATLAB
figs_autoarrange = true;  
lg = 0.01;  % 0.01 % logarithm offset for spectrum display
% KO Kx-Ky cross section for Kz:
% KO_z_crossect = 'ODT';

if length(unique(size(n_obj))) > 1
    error('`n_obj` must be a cube with equal dimensions.')
elseif ~isempty(n_obj) && size(n_obj,1)~=Nx
    error('`n_obj` must be of size Nx.');
end
if length(unique(size(x_tv))) > 1
    error('`mask_obj` must be a cube with equal dimensions.')
elseif ~isempty(x_tv) && size(x_tv,1)~=Nx
	error('`mask_obj` must be of size Nx.')
end
switch geometry
case 'fixed'
    if ~isempty(thetay); error(['Geometry `fixed` is not allowed in `thetay` ' ...
								'illumination scenario (`thetay` not empty).']); end
	if Ramp; warning('Using experimental ramp2D-filter for LAT (conical illumination).'); end
case 'facing'
	if isempty(thetay); warning(['Using arbitrary illumination scenario (kx,ky) ' ...
								'with geometry `facing` (not LAT).']);
						warning('<< Press any key to confirm >>'); pause();
	else; warning(['Assuming full-angle configuration around Y axis. ' ...
					  'Make sure that projections are oriented properly.']);
	end
otherwise
    error('Wrong option `geometry`.')
end
if size(Fpmask,3)==1
	Fpmask = repmat(Fpmask,[1 1 nproj]);
elseif sum(size(Fpmask(:,:,1))~=size(SINOamp(:,:,1)))
	error('Size of Fpmask does not match the sinogram.')
end

%%

kn = n_imm./sino_params(3,:); % length of wavevector

kxp        = sino_params(1,:);
kyp        = sino_params(2,:);
lambda_all = sino_params(3,:);
NA         = sino_params(4,:);
ind_oct    = find(sino_params(5,:)==2)';
ind_odt    = find(sino_params(5,:)==1)';

kxy = sqrt(kxp.^2 + kyp.^2);
fprintf('Max zenith angle: %0.2f\n', max(asind(kxy./kn)));
fprintf('Min zenith angle: %0.2f\n', min(asind(kxy./kn)));

% z-shift in K space due to xy shifts
% shifts of the center of Ewald spheres for all projections
% based on kxp and kyp, that are the locations of peaks within aperture
kz = real(sqrt( kn.^2 - kxp.^2 - kyp.^2 ));
if ~all(kz>=0); error('Some (kxp,kyp) describe evanescent waves.'); end

N_projection = round(N_projection/2)*2; 
if N_projection<Nx
	error('Padded projection size must be equal or larger than original size: NP>=Nx.');
end
dkP = 1/(dx_projection*N_projection); % frequency sample (padded projection)

% K space size
OTF_diameter=[];
if strcmp(geometry,'fixed')
    if strcmp(Kspace_padding,'optimal')
        k=0;
        for jj = 1:numel(ind_oct)
            ii = ind_oct(jj);
            k = k+1;
            % if xy-shift doesn't result in arc being outside of Ewald sphere, then only z-shift increases OTF-diameter
            if ( NA(ii)/lambda_all(ii)+sqrt(kxp(ii)^2+kyp(ii)^2) - kn(ii))<0
                % OTF_diameter = 2*( (Ewald sphere radius)+(z-shift of arc) )
                OTF_diameter(k) = 2*( kn(ii) + sqrt(kn(ii)^2-kxp(ii)^2-kyp(ii)^2) );
            % otherwise the OTF-diameter could be increased due to xy-shift or z-shift of arc
            else
                % check what is greater: (xy-extent of arc) + (xy-shift of arc) - (Ewald sphere diameter) OR (z-shift of arc)
                % in other words: which 'surplus' over Ewald sphere diameter is greater: in xy or z direction
                OTF_diameter(k) = 2*( kn(ii) + max(NA(ii)/lambda_all(ii)+sqrt(kxp(ii)^2+kyp(ii)^2)-kn(ii), sqrt(kn(ii)^2-kxp(ii)^2-kyp(ii)^2)) );
            end
        end
        if numel(ind_odt)
            OTF_diameter(end+1) = 2*(max(NA(ind_odt))/min(lambda_all(ind_odt)) + max( sqrt( kxp(ind_odt).^2+kyp(ind_odt).^2))); % 2*(NA/lambda + max existing xy shift)
        end
        OTF_diameter = max(OTF_diameter);

    elseif strcmp(Kspace_padding, 'maxaperture')
        if numel(ind_oct)
            OTF_diameter = 4*max(kn(:));
        else
            OTF_diameter = 4*(max(NA(ind_odt))/min(lambda_all(:))); % max aperture
        end
    end
end

if strcmp(geometry,'facing')  % from pole to equator of Ewald sphere limited by NA, no OCT here!
    OTF_diameter = 2*sqrt((max(kn)-sqrt(max(kn)^2-(max(NA(ind_odt))/min(lambda_all))^2))^2+(max(NA(ind_odt))/min(lambda_all))^2);
end

% N_Kspace_xy = min(max_matrix_size,round(OTF_diameter/dkP /2)*2);
N_Kspace_xy = round(OTF_diameter/dkP /2)*2;
% make it divisible by 12...48? (common for 12 and 16 (max crop 0.25 and max dwnsample 4=16))
N_Kspace_xy = round(N_Kspace_xy/48)*48;

N_Kspace_z_upsampled = Kspace_oversampling_z*N_Kspace_xy;
N_Kspace_z_upsampled_cropped = N_Kspace_z_upsampled*limit_resolution_z/ROI_crop_z;

if N_Kspace_xy==0; error('N_Kspace_xy_padded==0. Wrong units?'); end

% % Scale ROI_crop_z. If Kspace_oversampling_z is e.g. 2 and ROI_crop_z=2, then we dont want
% % the size of Kspace to go back to nominal, but to be actually 2x
% % downsampled compared to nominal:
% ROI_crop_z = ROI_crop_z*Kspace_oversampling_z;

dkPz = dkP/(N_Kspace_z_upsampled/N_Kspace_xy)*ROI_crop_z;

% projection spectrum plane
if strcmp(geometry,'fixed')
	NF = N_Kspace_xy; % must contain shadows of all shifted Ewald spheres
    % additionally upsample projection spectra for inpterpFp:
	NF = NF*q; N_projection = N_projection*q; dkP = dkP/q;
elseif strcmp(geometry,'facing')
	NF = N_projection; % original projection band (upsampled for NP>Nx)
end
ks = ((1:NF)-floor(NF/2)-1)*dkP; % zero as in fftshift(fft(x))
[kyy, kxx] = meshgrid(ks); % x-vertical

% K-space formation from projections stored in SINO arrays
disp(['Assembling K space from ', num2str(nproj), ' projection spectra'])
% K-space (FT of scattering potential):
KO = single(zeros(N_Kspace_xy,N_Kspace_xy,N_Kspace_z_upsampled_cropped)); 
EW = single(zeros(size(KO))); % voxel hit accumulator


if all(sino_params(5,:)==2)      
    Kz_slice = 50;
else
    Kz_slice = N_Kspace_z_upsampled_cropped/2+1;
end

% if nproj<10; jump = 1; 
% elseif nproj<50; jump = 5; 
% else jump=10; 
% end

for j=1:nproj  % main loop for K-space generation
    % Ewald sphere mesh centered around (0,0)
    rf = NA./lambda_all/dkP;
    
    NF0 = 2*round(rf(j)); 
    ks0 = ((1:NF0)-floor(NF0/2)-1)*dkP; 
    [kyy0, kxx0] = meshgrid(ks0);  % x-vertical

    kzz0 = real(sqrt(kn(j)^2 - kxx0.^2 - kyy0.^2)); 

%     if (mod(j,jump)==0 || j==nproj); disp([num2str(j) '/' num2str(nproj)]); end
    
    % define projection field
	amp = SINOamp(1:Nx,1:Ny,j); amp(amp==0 | isnan(amp)) = 1; % amp=0 not acceptable for Rytov
	ph = SINOph(1:Nx,1:Ny,j); ph(isnan(ph)) = 0;
    switch Approx
        case 'Born'
            Up = single(amp.*exp(1i*ph)-1);
        case 'Rytov'
            Up = single(log(amp)+1i*ph);     % log(amp * exp(1i*ph)) = log(amp)+1i*ph
        otherwise
            error('Wrong option `Approx`.');
    end
    
	% projection padding (does not matter if spectrum is wrapped, will get upsampled)
	if strcmp(Approx,'Born') % always bright-field, so Tukey before zeropadding
        Up = tukeywin(Nx, tukr)*tukeywin(Ny, tukr).' .* Up; % tukSa
    end
	Up = padarray(Up, [(N_projection-Nx)/2 (N_projection-Ny)/2], 0);% abs(Up) should be >=0, as SINOph; angle(Up) should oscilate around 0
	Fpm = round(imresize(single(Fpmask(:,:,j)),[N_projection N_projection]));
				
	% projection spectrum (NA filtering)
	if objshift
    	Fp = fftshift(fft2(ifftshift(Up)))*dx_projection^2;
    else
        Fp = fftshift(fft2(Up))*dx_projection^2;
    end
    
	if isempty(thetay) % LAT
		if N_projection<NF && N_projection>=2*NA(j)/lambda_all(j)/dkP-2% zip=1 (2*NA < projection band < 4*NA); rf = NA/lambda/dkP
			Fp = circshift(Fp, round([kxp(j) kyp(j)]/dkP)); % center NA (instead of energy peak) to unwrap (zip=1)
			Fp = subarray_yx(Fp, [NF NF], [], 0); % crop/pad to match kxx,kyy,kzz (X-Y slice of KO)
			Fp = circshift(Fp, round([-kxp(j) -kyp(j)]/dkP)); % shift back perfectly (peak was already centered with subpixel precision)
            % if~isempty(Fpmask0);warning('Fpmask ignored for compressed projection (zip=1).');end
            % Fpm = true(size(Fp));% Fpmask ill-defined for zip=1
            Fpm = circshift(Fpm, round([kxp(j) kyp(j)]/dkP)); % center NA to unwrap (zip=1)
			Fpm = subarray_yx(Fpm, [NF NF], [], 0); % crop/pad to match kxx,kyy,kzz
			Fpm = circshift(Fpm, round([-kxp(j) -kyp(j)]/dkP)); % shift back perfectly
            
		else % zip=0 or projection band <2*NA (narrow-band)
			Fp = subarray_yx(Fp, [NF NF], [], 0); % crop/pad to match kxx,kyy,kzz (X-Y slice of KO)
            Fpm = subarray_yx(Fpm,size(Fp), [], 1);
		end
        % Fpm = subarray_yx(Fpm,size(Fp), [], 1);
    end		
	
    % show spectrum of projection with dimmed out frequencies aroud NA
    if contains(plots,'DD')
        if isempty(thetay); offset = [-kxp(j) -kyp(j)]/dkP; else; offset = [0 0]; end
        Fp_=Fp.*(1-.95*(1-ellipse_yx(size(Fp),[rf(j) rf(j)],offset)));
        figure(10);imagesc(log10(lg+abs(Fp_.')));axis image;%colorbar
        xlabel('Kx'); ylabel('Ky'); title('Fp');colormap jet;drawnow
        % waitforbuttonpress
    end
     
	% Ewald insertion coords in K;
	% (Kx,Ky,Kz)=(0,0,0) corresponds to the power peak in Fp at (0,0)
	if strcmp(geometry,'fixed')
		if ~strcmp(interpFp,'none') % Fp interpolated onto centered Ewald mesh
            % Move the coordinate system associated with the Ewald sphere 
            % according to the position where this Ewald sphere should be
            % placed in the K-space.
            % This is required because in this configuration the OTF in the
            % Fp is centered
            Kx = kxx0 - kxp(j); 
			Ky = kyy0 - kyp(j);
			kzz = kzz0;
			Kz = kzz0 - kz(j);
            %------------------------------------------------------------------------------------------------
            if any(ind_oct==j)
                Kz = -(kzz0 + kz(j));  % <-- inverted arc as in reflection ODT
            end
            %------------------------------------------------------------------------------------------------
		else % Ewald sampled by Fp(NA) mesh
			Kx = kxx; % kxx==kxx0 due to q=1 for interpFp='none'
			Ky = kyy;
			kzz = real(sqrt(kn(j)^2-(kxx + kxp(j)).^2-(kyy + kyp(j)).^2 )); % irregular Ewald mesh above Fp(NA)
            Kz = kzz - kz(j);  % kzz - kz(j): shift of -kz to 0 (coord origin); Kz = kzz: no shift;
            %------------------------------------------------------------------------------------------------
            if any(ind_oct==j)
                Kz = -(kzz + kz(j)); % <-- inverted arc, shift of -3*kz_oct to the proper position
            end
        end
    
	elseif strcmp(geometry,'facing')
		Kx = kxx0;
		Ky = kyy0;
        Kz = kzz0 - kn(j);
		
		kzz = kzz0;

		if ~isempty(thetay) % scenario 'thetay'
			% rotate the c.s. with left-hand rule (around Ky), so that the data is rotated
			% by theta_x with right-hand rule (around Y):
			Kxr = Kz*sin(thetay(j)) + Kx*cos(thetay(j));
			Kzr = Kz*cos(thetay(j)) - Kx*sin(thetay(j));
			Kx = Kxr;
			Kz = Kzr;
		else % arbitrary scenario
			fid = atan2d(kyp(j), kxp(j));
			%thd = real(asind( min( 1, sqrt(kxp(j)^2+kyp(j)^2)/(n_imm/lambda) ) ));
            thd = real(asind( min( 1, sqrt(kxp(j)^2+kyp(j)^2) / kn(j) ) ));
            
			if abs(thd-pi/2)<1e-12; thd=pi/2-1e-12; end % avoid singularity
			RM = SpinCalc('EA123toDCM',-[0 0 fid],1e-6,false) * ...
			     SpinCalc('EA123toDCM',-[0 thd 0],1e-6,false); % theta-POLE
			K_XYZ = RM * [Kx(:).'; Ky(:).'; Kz(:).'];
			Kx = K_XYZ(1,:).';
			Ky = K_XYZ(2,:).';
			Kz = K_XYZ(3,:).';
		end
	end
	
    % contribution to K space
    if Ramp
		if strcmp(geometry,'fixed') % conical band-limited filter for conical illumination scenario
			RF = ramp_filter(ks,2);
		else % band-limited ramp filter (Marone & Stampanoni, 2012)
			RF = ramp_filter(ks,1); % column
            % % RF = RF .* [0; parzenwin(length(RF)-1)]; % supress high-freq noise (max at ks=0)
            % RF = RF .* parzenwin(length(RF)); % assymetric, but better?
			RF = repmat(reshape(RF,[],1),[1 NF]); % X-vertical, rotation aroud Y
		end
		Fp = Fp.*RF;
    end
    
	% shifting Fp content under symmetrical Ewald mesh
	if strcmp(geometry,'fixed')
        if strcmp(interpFp,'none')
            % Fp=Fp;
        elseif strcmp(interpFp,'linc')
            if objshift; Up = fftshift(ifft2(ifftshift(Fp)))/dx_projection^2; FpmFT = fftshift(ifft2(ifftshift(Fpm)));
            else;		 Up = ifft2(ifftshift(Fp))/dx_projection^2; FpmFT=ifft2(ifftshift(Fpm)); end
            Up = lincarrier2d(Up, [kxp(j) kyp(j)]/dkP);
            FpmFT = lincarrier2d(FpmFT, [kxp(j) kyp(j)]/dkP);
            % prevent "wiskers" in spectrum (experimental data: visible side effects!)
            Sa = abs(Up);
            av = 0; % 0! because phase fringes are discontinuous on edges
            Sa = tukeywin(NF, tukr)*tukeywin(NF, tukr).' .* (Sa - av) + av;
            Up = Sa.*exp(1i*angle(Up));
            if objshift; Fp = fftshift(fft2(ifftshift(Up)))*dx_projection^2; Fpm = fftshift(fft2(ifftshift(FpmFT)));
            else;		 Fp = fftshift(fft2(Up))*dx_projection^2;	Fpm = fftshift(fft2(FpmFT)); end
            Fp = subarray_yx(Fp, size(kzz));
            Fpm = subarray_yx(abs(Fpm)>0.5,size(Fp),[],1);
        else% (nearest, linear, cubic, spline)
            % interpolate fragment of Fp to find projection spectrum values directly
            % below the symmetrical Ewald mesh:
            Fp = interp2(kyy,kxx, Fp, Ky,Kx, interpFp); % high-quality interpolation required
        end
	end
	if strcmp(geometry,'facing')
		Fp = subarray_yx(Fp, size(kzz)); % match size with symmetric Ewald mesh
		Fpm = subarray_yx(Fpm,size(Fp),[],0);
	end

    % generate a circle (Fpm) that represents OTF of a projection 
    % (-> limited by NA of microscope objective)
    if strcmp(geometry,'fixed') && strcmp(interpFp,'none') 
        Fpm = Fpm .* ellipse_yx(size(Fp), [rf(j) rf(j)], [-kxp(j) -kyp(j)]/dkP);
    else % 'facing' or 'fixed' with interpolated Fp
        Fpm = Fpm .* ellipse_yx(size(Fp), [rf(j) rf(j)]);
    end
    % multiply the Fourier spectrum of a projection (Fp) with a generated
    % OTF of a projection (Fpm)
    Fp = Fp .* Fpm;

    % show projection spectrum to be mapped on a sphere
    if contains(plots,'2')
        % ---------------------------------------------------------------------
        figure(11);imagesc(log10(lg+abs(Fp.')));axis image; %colorbar
        xlabel('Kx');ylabel('Ky');title('Fp  (to be mapped on a sphere)');colormap jet;drawnow  
        % ---------------------------------------------------------------------
    end
    % FFT2{proj(x,y)} -> FFT3{o(x,y,z)}(Ewald) FT of the scattering potential
    Fp = +1i.*(2*pi)*2*kzz.*Fp; % to be mapped on a sphere

	% discard evanescent waves and waves exceeding NA
	iiF = find(kzz>0 & Fpm~=0);
    % floating point subscripts addressing N_Kspace_xy_padded x N_Kspace_xy_padded x N_Kspace_xy_paddedz array
    Kx = Kx/dkP/q+ceil((N_Kspace_xy+1)/2); % N_Kspace_xy_padded / 2 + 1 corresponds to position of the DC term
    Ky = Ky/dkP/q+ceil((N_Kspace_xy+1)/2);
    Kz = Kz/dkPz+ceil((N_Kspace_z_upsampled_cropped+1)/2);
    % discard frequencies exceeding N_Kspace_xy_padded x N_Kspace_xy_padded x N_Kspace_z_padded_upsampled 
    % array subscripts
    iiF( Kx(iiF)<1 | Kx(iiF)>N_Kspace_xy ) = [];
    iiF( Ky(iiF)<1 | Ky(iiF)>N_Kspace_xy ) = [];
    iiF( Kz(iiF)<1 | Kz(iiF)>N_Kspace_z_upsampled_cropped ) = [];
		
    % rounded subscripts (nearest-neighbour mapping)
    Kx = round(Kx(iiF));
    Ky = round(Ky(iiF));
    Kz = round(Kz(iiF));    
    % linear indices of K space points covered by projection data
    % (wraps spectra in XY plane with single-voxel jumps in Z, if N_Kspace_xy_padded is too small)
    iiK = (Kz-1)*N_Kspace_xy^2 + (Ky-1)*N_Kspace_xy + Kx; 
    % iiK = sub2ind([N_Kspace_xy_padded N_Kspace_xy_padded N_Kspace_xy_padded], Kx, Ky, Kz); % returns an error if N_Kspace_xy_padded is too small

    KO(iiK) = KO(iiK) + Fp(iiF); % mapping Fp onto Ewald sphere or plane!
    EW(iiK) = EW(iiK) + 1; % voxel hit counter ++

%     if contains(KO_z_crossect, 'ODT')
%         Kz_slice = floor(size(KO, 3)/2)+1;
%     elseif contains(KO_z_crossect, 'OCT') 
%         Kz_slice = 50;
%     end
    
    % show the process of filling K-space with projections
    if contains(plots,'D')
        fig20 = figure(20);
        if ~verLessThan('matlab', '9.5'); sgtitle( sprintf('KO (%d/%d)',j,nproj) ); end
        %
        subplot(121); 
        imagesc(log10(lg+abs(squeeze(KO(:,end/2+1,:)).'))); axis square;
        xlabel('Kx'); ylabel('Kz'); title('KO  Kx-Kz plane');
        %
        subplot(122); 
        imagesc(log10(lg+abs(squeeze(KO(:,:,Kz_slice)).'))); axis square;
        xlabel('Kx'); ylabel('Ky'); title('KO  Kx-Ky plane');
        colormap jet
        drawnow
    end
     
    h = findobj('type','figure');
    figs_number = length(h);
    if figs_number==1 && exist('fig20', 'var') 
        fig20.set('Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
        clear fig20;
    elseif figs_number>1 
        if figs_autoarrange; autoArrangeFigures(); end
    end
    clear h figs_number
    
end

% clear unused variables
clear Fpmask Fpmask0

if ~strcmp(interpFp,'none')
	N_projection=N_projection/q;dkP=dkP*q;% back to actual values
end
if ~Ramp	
	KO(EW>0) = KO(EW>0) ./ single(EW(EW>0));% average out multiple hits in voxels
end

EW = (EW>0); % voxels with data

clear SINOph SINOamp

%make the spectrum hermitian
% KO = makehermitian(KO);

% O(Kx,Ky,Kz) --> n(x,y,z)
dxo_xy = dx_projection*N_projection/N_Kspace_xy;% sample size in reconstruction space (NP/N_Kspace_xy_padded: how many more dk samples in K)
dxo_z = dxo_xy/limit_resolution_z;

if ~isempty(x_tv)
	disp('Signal initialization from total-variation reconstruction ...')	
	if any(isnan(x_tv)); error('NaNs detected in `init_obj`.'); end

    x_tv = padarray(x_tv, [1 1 1]*(N_projection-Nx)/2, 0); % pad same way as projections
    x_tv = ndresize(x_tv, [N_Kspace_xy N_Kspace_xy N_Kspace_xy*limit_resolution_z], 'linear'); % match spectrum padding
    if Kspace_oversampling_z>1
        x_tv = padarray(x_tv, [0 0 (N_Kspace_xy*limit_resolution_z*Kspace_oversampling_z-N_Kspace_xy*limit_resolution_z)/2], 0); % match padding from spectrum oversampling in Z
    elseif ROI_crop_z >1
        x_tv = ndcrop(x_tv, [size(x_tv,1) size(x_tv,2) size(KO,3)]);
    end
    %% Create object support from x_tv
    object_support = fcreatemask3d(abs(x_tv),obj_type);
    object_support = imfill(object_support,'holes');

    n_rec = x_tv+n_imm;		
    if objshift; n_rec = ifftshift(n_rec); object_support = ifftshift(object_support); end% mimic shape from ifftn
else
    object_support = [];
	disp('Direct inversion ...')
       
    n_rec = ifftn(ifftshift(KO))./(dxo_xy^2*dxo_z);
	if ~objshift; n_rec = fftshift(n_rec,3); end% remove z-shift
    % Born, Rytov
    n_rec = n_imm*sqrt(1-n_rec/(2*pi*mean(kn))^2); % o(x,y,z) -> n(x,y,z)
end

% show fully filled K-space, Direct Inversion result and object support (if applies)
if contains(plots,'3')
    if objshift; n_rec_0 = fftshift(n_rec); 
    else		 n_rec_0 = n_rec; end
    figure(30); % set(figure(30),'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
    if ~verLessThan('matlab', '9.5'); sgtitle('KO - before GP'); end
    subplot(221);
    imagesc(log10(lg+abs(squeeze(KO(:,end/2+1,:)).'))); axis square;  
                  xlabel('Kx'); ylabel('Kz'); title('KO  Kx-Kz plane');
    subplot(223);
    %imagesc(log10(lg+abs(squeeze(KO(:,:,end/2+1)).'))); axis image;  
    imagesc(log10(lg+abs(squeeze(KO(:,:,Kz_slice)).'))); axis image;  
                  xlabel('Kx'); ylabel('Ky'); title('KO  Kx-Ky plane');
    
    n_rec_0_temp = ndcrop(n_rec_0,[N_Kspace_xy,N_Kspace_xy,N_Kspace_xy]);
    n_rec_0_temp_xz = real(squeeze(n_rec_0_temp(:,end/2+1,:)).');
    n_rec_0_temp_xy = real(squeeze(n_rec_0_temp(:,:,end/2+1)).');
    
    if ~isempty(object_support)
        if objshift
            object_support_temp = ndcrop(fftshift(object_support),[N_Kspace_xy,N_Kspace_xy,N_Kspace_xy]);
        else
            object_support_temp = ndcrop(object_support,[N_Kspace_xy,N_Kspace_xy,N_Kspace_xy]);
        end
        object_support_xz = squeeze(object_support_temp(:,end/2+1,:)).';
        object_support_xz_contour = abs(del2(object_support_xz>0.5));
        n_rec_0_temp_xz(object_support_xz_contour>0) = max(n_rec_0_temp_xz(:));
        
        object_support_xy = squeeze(object_support_temp(:,:,end/2+1)).';
        object_support_xy_contour = abs(del2(object_support_xy>0.5));
        n_rec_0_temp_xy(object_support_xy_contour>0) = max(n_rec_0_temp_xy(:));
    end
    
    subplot(222); imagesc( n_rec_0_temp_xz ); axis image; colormap jet   
                      xlabel('x'); ylabel('z'); 
    if ~isempty(object_support)
        title('Direct Inversion with object support contour (XZ)');
    else
        title('Direct Inversion (XZ)'); 
    end
    
    subplot(224); imagesc( n_rec_0_temp_xy ); axis image; colormap jet
                      xlabel('x'); ylabel('y');
    if ~isempty(object_support)
        title('Direct Inversion with object support contour (XY)');
    else
        title('Direct Inversion (XY)'); 
    end

    drawnow
    clear n_rec_0_temp n_rec_0_temp_xz n_rec_0_temp_xy object_support_xz ...
        object_support_xz_contour object_support_xy object_support_xy_contour
end

% GP settings
RMAEtab = single(NaN(1,nGPi+1));
RRMSEtab = single(NaN(1,nGPi+1));
RMADtab = single(NaN(1,nGPi+1));
RRMSDtab = single(NaN(1,nGPi+1));
% define object block surrounded by borders from projection padding
block = subarray_yx(true(1,round(N_Kspace_xy/projection_padding_xy/2)*2),[1,N_Kspace_xy]);
blockz = subarray_yx(true(1,round(N_Kspace_z_upsampled_cropped/projection_padding_xy/2)*2),[1,N_Kspace_z_upsampled_cropped]);
qmask = true(nnz(block),nnz(block),nnz(blockz));
if objshift; block = ifftshift(block);
			 blockz = ifftshift(blockz); end
if ~isempty(n_obj)
	n_rec_0 = n_rec;
	% object padding to match projection padding
	n_obj_0 = padarray(n_obj, [1 1 1]*(N_projection-Nx)/2, n_imm); % pad same way as projections (n_imm!)
    % match spectrum padding:
	n_obj_0 = ndresize(n_obj_0, [N_Kspace_xy N_Kspace_xy N_Kspace_xy], 'linear');
    % match padding from spectrum oversampling in Z:
	n_obj_0 = padarray(n_obj_0, [0 0 (N_Kspace_z_upsampled-N_Kspace_xy)/2], n_imm); 
	% mask for QI evaluation
	if do_qmask; qmask = imfill((abs(n_obj_0-n_imm)>0),'holes');
	else		 qmask = true(size(n_obj_0));end
	% mimic current shape of n_rec
	if objshift; n_obj_0 = ifftshift(n_obj_0);
				 qmask = ifftshift(qmask); end
	n_rec_0 = n_rec_0(block,block,blockz);
	n_obj_0 = n_obj_0(block,block,blockz);
	qmask = qmask(block,block,blockz);

	RMAEtab(1) = rpnorm(real(n_rec_0),real(n_obj_0),1);
	RRMSEtab(1) = rpnorm(real(n_rec_0),real(n_obj_0),2);
	RMADtab(1) = NaN;
	RRMSDtab(1) = NaN;
end
if nGPi==0
	KOi = KO;% should not duplicate memory
end

% replace 3D array with lists of occupied indices and values:
EWi = find(EW>0); EKO = KO(EWi); clear KO;

if relaxLaR>0
	cover_out = subarray_3d(single(ones(round([N_Kspace_xy/projection_padding_xy, N_Kspace_xy /... 
                                          projection_padding_xy, N_Kspace_xy/projection_padding_xy]/2)*2)), ...
                                          [N_Kspace_xy N_Kspace_xy N_Kspace_z_upsampled]);% LaR(best:Shepp80?)
	if objshift; cover_out = ifftshift(cover_out); end% mimic current shape of n_rec
	Pi = find(cover_out~=1);% (blurred) padding region
	Pc = cover_out(Pi);% padding mask
	betaLaR = 1;
	clear cover_out
end

if nGPi>0;disp('GP iterations ...');end
if nGPi<20; jump = 1;
else;       jump = 5; end

for gpi = 1:nGPi
    
    if mod(gpi,jump)==0; disp([num2str(gpi) '/' num2str(nGPi)]); end
    
	n_rec_prev_ev = ndresize(n_rec(block,block,blockz), [1 1 1]*min(nnz(block),Nc), 'linear');
	% transparency constraint
	if do_TC
		n_rec = real(n_rec);
	end	
    % non-negativity constraint
	if do_NNC
		if do_NNC>0; id = (real(n_rec)<n_imm);% non-negativity
		else;		 id = (real(n_rec)>n_imm); end% non-positivity
        % n_rec(id) = n_imm + real(n_rec(id)-n_imm).*(1-cover_out(id)) + 1i*imag(n_rec(id)); % prev
		n_rec(id) = n_imm + 1i*imag(n_rec(id));
	end	
	
	% padding suppression (LaRoque)
	if relaxLaR>0 % (LaR,TUK-LaR)
       % n_rec = (n_rec-n_imm).*cover_out + ...
       %			  (zeros(size(n_rec))*betaGP + (n_rec-n_imm)*(1-betaGP)).*(1-cover_out) + n_imm;
		n_rec(Pi) = (n_rec(Pi)-n_imm).*((Pc-1)*betaLaR+1) + n_imm;
		betaLaR = relaxLaR*betaLaR;
	end	
	
    % spatial support
    if ~isempty(object_support)
        n_rec = (n_rec-n_imm).*(object_support*betaM+(1-betaM)) + n_imm;
        betaM = relaxM*betaM;
    end
	
    % O(Kx,Ky,Kz) <-- n(x,y,z)
    n_rec = (2*pi*mean(kn))^2*(1-n_rec.^2/n_imm^2);

	if ~objshift; n_rec = ifftshift(n_rec,3); end
    KOi = fftshift(fftn(n_rec)).*(dxo_xy^2*dxo_z);
    % replenishment with projection data
    % KOi = EW.*(KO*betaGP + KOi*(1-betaGP)) + (single(1)-EW).*KOi;% relax-inefficient
    % KOi(EWi) = EW(EWi).*(EKO*betaGP + KOi(EWi)*(1-betaGP)) + (single(1)-EW(EWi)).*KOi(EWi);% relax-efficient
    %
	KOi(EWi) = EW(EWi).*EKO + (single(1)-EW(EWi)).*KOi(EWi);% relax only delays convergence
    betaGP = relaxGP*betaGP;       
    
	% O(Kx,Ky,Kz) --> n(x,y,z)
    n_rec = ifftn(ifftshift(KOi))./(dxo_xy^2*dxo_z);
	if ~objshift; n_rec = fftshift(n_rec,3); end
    n_rec = n_imm*sqrt(1-n_rec/(2*pi*mean(kn))^2); % o(x,y,z) -> n(x,y,z)
	
	n_rec_ev = ndresize(n_rec(block,block,blockz), [1 1 1]*min(nnz(block),Nc), 'linear');
    % RMADtab(gpi+1) =  norm(real(n_rec_ev(:)-n_rec_prev_ev(:)),1)/norm(real(n_rec_prev_ev(:)),1);
    % RRMSDtab(gpi+1) = norm(real(n_rec_ev(:)-n_rec_prev_ev(:)),2)/norm(real(n_rec_prev_ev(:)),2);
    RMADtab(gpi+1) =  rpnorm(real(n_rec_ev),real(n_rec_prev_ev),1);
    RRMSDtab(gpi+1) = rpnorm(real(n_rec_ev),real(n_rec_prev_ev),2);
	
    % show reconstruction (XZ) progress in GP and 1D plot of relative progress of GP
    if contains(plots,'4')
        figure(35)
        n_rec_0 = n_rec(block,block,blockz);
        if objshift; rec_xz = fftshift(real(squeeze((n_rec_0(:,1,:)-n_imm).*qmask(:,1,:)+n_imm).'));
        else;		 rec_xz = real(squeeze((n_rec_0(:,end/2+1,:)-n_imm).*qmask(:,end/2+1,:)+n_imm).'); end
        if isempty(n_obj);subplot(211);else;subplot(221);end
        imagesc(rec_xz); axis image; %colorbar % [1.55 1.6]
        xlabel('x'); ylabel('z')
        title(sprintf('(%d/%d)',gpi,nGPi))
        if isempty(n_obj);subplot(212);else;subplot(223);end
        plot(0:nGPi, RMADtab, 'k-x');hold on
        plot(0:nGPi, RRMSDtab, '-+');hold off
        grid on
        % ylim([0 5e-3]);xlim([0 length(RMADtab)])
        legend({'RMAD','RRMSD'},'Location','northeast')
        % title(sprintf('delta=%g > epsilon=%g', abs(RRMSDtab(gpi+1)-RRMSDtab(gpi)),epsi))
        title(sprintf('delta=%g > epsilon=%g', abs(RMADtab(gpi+1)-RMADtab(gpi)),epsi))
        colormap jet
        drawnow
        if gpi==1;set(gcf,'Position',[0 0 1920/2 1080]);end
    end	
	
	% stopping conditions
	% (increment-based): RMADtab
	if ~isnan(epsi) && gpi>2 && abs(RMADtab(gpi+1)-RMADtab(gpi)) < epsi% absolute slope
		disp([num2str(gpi) '/' num2str(nGPi)]);
		fprintf('Stop on RMAD saturation (epsilon=%g).\n',epsi)
		nGPi = gpi;
		break;
    end
    
    if gpi==nGPi; disp([num2str(gpi) '/' num2str(nGPi)]); break; end
    
    
%     if objshift; RECON = fftshift(n_rec); end
%     RECON = permute(real(RECON), [2 1 3]); % (y,x,z)
%     % N_Kspace_z_padded_upsampled cropped to N_Kspace_xy_padded / alfa (cropping the Z axis to match sizes in X and Y directions):
%     RECON = ndcrop(RECON, round( [N_Kspace_xy_padded N_Kspace_xy_padded N_Kspace_xy_padded] / projection_padding_xy /2)*2);
% 
%     img1 = squeeze(RECON(:,size(RECON,2)/2,:));
%     img2 = squeeze(RECON(:,:,size(RECON,3)/2));
%     figure(78)
%     subplot(1,2,1), imagesc(rot90(img1),[1.33 1.37]), axis image
%     subplot(1,2,2), imagesc(      img2 ,[1.33 1.37]),  axis image
%     a = axes;
%     t = title(num2str(gpi,'%03.f'));
%     a.Visible = 'off';
%     t.Visible = 'on';
%     drawnow
%     saveas(gcf,['/run/media/wojtek/Data2/SynologyDrive/Data/PROC_03_hacat_03/gpfos/','recon',num2str(gpi,'%03.f'),'.png'])
    
end

% apply constraints to final result
% if nGPi>0 && do_constraint_rec
if do_constraint_rec
	if do_TC && nGPi>0 % leave imag for DI
		n_rec = real(n_rec);
	end
	if do_NNC && nGPi>0
		if do_NNC>0; id = (real(n_rec)<n_imm);
		else;		 id = (real(n_rec)>n_imm); end % "non-positivity" contraint
		n_rec(id) = n_imm + 1i*imag(n_rec(id));
	end
	if ~isempty(object_support)
 		n_rec = (n_rec-n_imm).*object_support + n_imm;
	end
	if relaxLaR>0 % (LaR,TUK-LaR)
        % n_rec = (n_rec-n_imm).*cover_out + ... 
        %			    (zeros(size(n_rec))*betaGP + (n_rec-n_imm)*(1-betaGP)).*(1-cover_out) + n_imm;
		n_rec(Pi) = (n_rec(Pi)-n_imm).*((Pc-1)*betaLaR+1) + n_imm;
	end
end


% correct for shifts
if objshift; n_rec = fftshift(n_rec); end
% ramp-related normalization factor (depends on nproj)
if Ramp 
	n_rec = n_rec - n_imm;
	n_rec = n_rec*pi/(2*nproj)/dx_projection /dkP;
	n_rec = n_rec + n_imm;
end

% free some memory and recreate KO
if exist('n_obj','var'); clear n_obj; end
if exist('mask_obj','var'); clear mask_obj; end
if exist('n_rec_0','var'); clear n_rec_0; end
KO = single(zeros(size(KOi)));
KO(EWi) = EKO;

% RECONSTRUCTION (y,x,z)
RECON = permute(real(n_rec), [2 1 3]); % (y,x,z)
% N_Kspace_z_padded_upsampled cropped to N_Kspace_xy_padded / alfa (cropping the Z axis to match sizes in X and Y directions):
RECON = ndcrop(RECON, [N_Kspace_xy/projection_padding_xy N_Kspace_xy/projection_padding_xy N_Kspace_z_upsampled_cropped/Kspace_oversampling_z]);
