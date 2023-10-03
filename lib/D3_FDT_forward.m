% -------------------------------------------------------------------------
% Forward projector for 3D Diffraction Tomography
% Copyright: 2015-2018, Piotr L. Makowski, p.makowski@mchtr.pw.edu.pl
% -------------------------------------------------------------------------
%
%   * projection FTs mapped from 3D object Fourier space (K space) as:
%		- Ewald spheres ('Born' or 'Rytov' projections: first order scattering approx)
%		- planes ('PhaseRay' projections: integrated phase)
%   * projection geometry:
%		- 'facing': detector perpendicular to illumination
%		- 'fixed': detector always horizontal
%	* illumination scenario:
%		- full-angle with object rotating around Y (thetay)
%		- limited-angle with arbitrary illumination around Z (kxp,kyp)
%   * linear/cubic interpolation of K space before mapping to projection spectra
%	* projection aliasing prevented by object padding
%
% Use 'D3_FDT' for object reconstruction.
%
% _______________________________________________________________________
% TODO: geometry=='fixed' with Ewald=0 (?)
% _______________________________________________________________________
% 29.08.2019
% - add zpc<1 option (object pre-padding opc=1/zpc)
% - prevent object downsampling
% 09.08.2019
% - force obj_sign to be real to prevent complex SINOph
% - correct peojection energy to unity
% - add option sinterp='map' to match n-n mapping from solver
% 13.08.2018
% - add argument Nxp0, make savemem automatic for proj bands <2*NA/lambda
% 24.07.2018
% - do not use subbaray_yx so that the smudges are not cropped assymetrically in
%   case of unchanged bandwidth (zip=0)
% - fix interp3 failure in Octave (interp3 reads only real part of the input array)
% 28.02.2018
% - fix bug with limited NA in case of geometry='fixed'
% - add internal params 'savemem' and 'Nxp0' for efficient computation of band-limited projections


function [SINOamp,SINOph, Nxp,dxp] = D3_FDT_forward(n_obj, mask_obj, kxp,kyp, thetay, ...
                           lambda,n_imm,NA,dxo, ...
						   geometry,Ewald,Approx, sinterp, ...
						   zip,Nxp0,zpc, ...
						   showplots)

%%%%%%%%%%%% INPUT
% n_obj         = 3D matrix (x,y,z) of the object refractive index
% mask_obj      - binary 3D array (x,y,z) for object masking before forward
%                 projection
% kxp, kyp      - limited-angle illumination scenario given by components
%                 of illuminating "wave vectors" |k|=n_imm/lambda,
%                 coresponding to projections in SINO(x,y,:);
%                 ignored if ~isempty(thetay)
% thetay        - full-angle illumination scenario given by array of
%                 projection angles measured around Y axis, coresponding to
%                 projections in SINO(x,y,:);
% lambda        - wavelength in vacuum
% n_imm         - refractive index of object immersion
% NA            - numerical aperture upper limit; =NaN gives NA=n_imm
% dxo           - sample size of n_obj
% zpc           - object padding coefficient to reduce projection aliasing
%				  (with Ewald=1) as well as improve spectrum interpolation
%				  (with Ewald=0,Approx='PhaseRay' and object tightly fitting the cube);
%				  projection dimension still corresponds to object cube dimension; 
%				  to obtain bigger projections, pass n_obj padded manually and set zpc=1;
%				  zpc<1 results in projection dimension reduced with respect to n_obj dimension
% geometry      - capture system configuration (detector orientation):
%                 'facing' = rotating object (projection angles around Y given in thetay)
%                 'fixed' = rotating illumination (projection directions given by kxp,kyp)
% Ewald         - true: map projectrion spectra on Ewald spheres instead of planes
%                 (located in K space in accordance to 'geometry' scheme)
% Approx        - 'Born', 'Rytov', 'PhaseRay' (projection phase == integrated phase)
% sinterp       - 'nearest', 'linear', 'cubic'(BEST), 'spline'
% zip			- (true): Compress projection band to one Ewald sphere (NA) to save sinogram 
%				  storage space. In case of geometry='fixed' projection spectra will be 
%				  wrapped within the spectral window without overlaps, which corresponds 
%				  to reversable downsampling in the signal domain. D3_FDT solver 
%				  automatically unwrapps projection spectra, therefore zip=1 has negligible
%				  impact on reconstruction quality.
% Nxp0			- (~=NaN): manually set projection spectrum array size for base 
%                 projection dimension (without zpc-related padding);
%				  the tails of the projection spectra are wrapped by default; if Nxp0 limits 
%				  projection bandwidth below 2*NA/lambda (savemem=true is triggered), wrapping is
%                 disabled and the tails are lost;the resulting sinogram should be reconstructed 
%                 with NK=NP to prvent allocation of 3D spectrum space for nonexisting tails
% showplots     - plots flikering annoyingly (debug mode)
%
%%%%%%%%%%%% OUTPUT
% SINOamp       - absolute value of projection field complex amplitudes (x,y,proj)
% SINOph        - unwrapped phase of projection field complex amplitudes (x,y,proj)
% Nx0           - projection size (square)
% dx0           - projection sample size


% -------------------------------------------------------------------------
% Apply fftshift projections before FFT to get smooth phase for Fourier mapping
objshift = 1;
% works with newer versions of MATLAB
figs_autoarrange = 1;
% logarithm offset for spectrum display
lg = 0.01;
% -------------------------------------------------------------------------
if length(unique(size(n_obj))) > 1
    error('Object matrix `n_obj` must be a cube with equal dimensions.')
end
if zpc<1
	warning('zpc<1 given: assuming pre-padding opc=1/zpc already present in the object.');
	opc = 1/zpc;% object pre-padding, projection frame dimension will be 1/zpc of object cube dimension, Nxp0 determines band
    zpc = 1;
else
	opc = 1;
end


switch geometry
case 'fixed'
	if ~Ewald; error('Geometry `fixed` not tested with Ewald=false (physical meaning?)'); end % TODO
    if ~isempty(thetay); error('geometry=`fixed` not compatible with non-empty `thetay`'); end
    % +1: shift projection spectra by [kxp,kyp]
    % -1: shift projection spectra by [-kxp,-kyp]
    kshift = 1; % +1
case 'facing'
    kshift = 0;
otherwise
    error('Wrong option `geometry`.')
end
% projection wave vector coords (rows)
kn = n_imm/lambda;
kxs = kshift*kxp(:).'; % rows
kys = kshift*kyp(:).';
kz = real(sqrt( kn^2 - kxs.^2 - kys.^2 ));
if ~all(kz>0)
	error('Some (kxp,kyp) describe evanescent waves. Wrong units?')
end
if isnan(NA)
	NA = n_imm;
elseif NA>n_imm
	error('NA must by <= n_imm')
end


% object
Nxo = size(n_obj,1); % original object
if mod(Nxo,2) 
	if ~isempty(mask_obj)
		if any(size(n_obj)-size(mask_obj))
			error('`mask_obj` does not match `n_obj` in size')
		end
		mask_obj = mask_obj(1:end-1, 1:end-1, 1:end-1);
	end
	n_obj = n_obj(1:end-1, 1:end-1, 1:end-1);
	Nxo = size(n_obj,1); % original object (even)
end
dko = 1/(dxo*Nxo); % original object frequency sample


% K space containing all shifted/rotated Ewald spheres
Bk0 = 2*NA/lambda; % Ewald bandwidth
if Nxp0<2+round(Bk0/dko /2)*2; savemem = true;% skip extracting irrelevant frequencies
else 						   savemem = false;end	
if strcmp(geometry,'fixed')
    k_max = max([kxs,kys]); % maximum Ewald sphere shift
    Bk = Bk0 + (1-savemem)*2*k_max;
elseif strcmp(geometry,'facing')
    Bk = Bk0*sqrt(2); % from pole to equator x 2
end
% spectrum padding to scenario bandwidth
% Resampling of the object is disabled in the case of existing pre-padding
% (opc) to prevent ambiguous downsizing problem when calculating Nx. It is
% assumed that the pre-padded object array has sufficient resolution to 
% contain complete set of Ewald spheres.
if opc==1
    NKo = 4 + round(Bk/dko /2)*2;% (4+ make sure the whole Ewald fits)
else
    NKo = Nxo;
end
NKo = max(Nxo,NKo);% prevent object downsampling
dxu = dxo*Nxo/NKo; % upsampled
% object padding for anti-aliasing
NKP = round(NKo*zpc /2)*2; 
dkP = 1/(dxu*NKP); % frequency sample (padded projection)


% projection spectrum plane matching K space dimensions (shadow)
NF = NKP;
ks = ((1:NF)-floor(NF/2)-1)*dkP; % zero as in fftshift(fft(x))
[kyy, kxx] = meshgrid(ks); % x-vertical
% coordinate transformation plan for extracted projection
% projection crop (remove anti-aliasing extension)
Nx = ceil(NKo/opc /2)*2 % Nx = ceil(NF*dkP/dko /2)*2
dx = 1/(dkP*NF); % =dxu if NF=NKP, not affected by crop
dk = 1/(dx*Nx); % ~dko due to Nx rounding [edit: ==dko!]
% target projection size and sample
if isnan(Nxp0)
	if zip
		Nxp = 2+ceil(round(Bk0/dk /2)/opc)*2;% (2+ make sure the whole Ewald fits)
		dxp = 1/(dk*Nxp);
	else
		Nxp = ceil(Nxo/opc /2)*2;% ~dko % original object bandwidth
%		Nxp = ceil(NKo/opc /2)*2;% optimal object bandwidth (assures NA circles will fit)
		dxp = 1/(dk*Nxp);% ~dxo
	end
else
	Nxp = Nxp0;
	dxp = 1/(dk*Nxp);
end


% (zpcFast: upsampling first, then padding)
obj_sign = sign(sum(real(n_obj(:))-n_imm));
n_obj = ndresize(n_obj, [NKo NKo NKo], 'linear'); ppad = (NKP-NKo)/2;
n_obj = padarray(n_obj, [ppad ppad ppad], n_imm);
% 	figure;imagesc(squeeze(real(n_obj(:,end/2,:))).');axis image % !!!!!!!!!!!!!!!!!!!!
if ~isempty(mask_obj)
	mask_obj = ndresize(mask_obj, [NKo NKo NKo], 'linear');
	mask_obj = padarray(mask_obj, [ppad ppad ppad], 0); % zeros!	
    if showplots
    % ---------------------------------------------------------------------
    figure(11)
    subplot(121);   imagesc(squeeze(mask_obj(:,end/2+1,:)).'); axis image
                    xlabel('x'); ylabel('z'); title('obj-mask padded to ZP1 & resized to ZPK')
    subplot(121);   imagesc(squeeze(mask_obj(end/2+1,:,:)).'); axis image
                    xlabel('y'); ylabel('z');
					colormap jet
	drawnow
    % ---------------------------------------------------------------------
    end
end


% object masking
if ~isempty(mask_obj)
    n_obj = (n_obj-n_imm).*mask_obj + n_imm;
end

if showplots
% -------------------------------------------------------------------------
figure(12)
subplot(222); imagesc(real(squeeze(n_obj(:,end/2+1,:)).')); axis image;
              xlabel('x'); ylabel('z')
              title('FDT forward: input object (masked)')
subplot(224); imagesc(real(squeeze(n_obj(:,:,end/2+1)).')); axis image;
              xlabel('x'); ylabel('y')
			  colormap jet
drawnow
% -------------------------------------------------------------------------
end

if strcmp(Approx,'PhaseRay')
    n_obj = n_obj - n_imm;
    n_obj = n_obj*(2*pi/lambda); % guessed factor (refractive index --> integrated phase)
else
    n_obj = (2*pi*kn)^2*(1-n_obj.^2/n_imm^2); % n(x,y,z) -> o(x,y,z) (scattering potential)
end
if objshift
    KO = fftshift(fftn(ifftshift( single(n_obj) )))*dxu^3; % O(x,y,z) = FFT3{o(x,y,z)}
else
    KO = fftshift(fftn(ifftshift( single(n_obj) ,3)))*dxu^3; % O(x,y,z) = FFT3{o(x,y,z)}
end
% 	KO(:,:,1:2) = 0;% @@@
% 	KO_xy = squeeze(KO(end/2+1,:,end:-1:1));% @@@
% 	KO_xy = squeeze(KO(:,:,end/2+1));% @@@

if showplots
% -------------------------------------------------------------------------
figure(12)
subplot(221); imagesc(log10(lg+abs(squeeze(KO(:,end/2+1,:)).'))); axis image;
              xlabel('Kx'); ylabel('Kz')
subplot(223); imagesc(log10(lg+abs(squeeze(KO(:,:,end/2+1)).'))); axis image;
              xlabel('Kx'); ylabel('Ky')
			  colormap jet
drawnow
% -------------------------------------------------------------------------
end


if ~isempty(thetay); nproj = length(thetay);
else                 nproj = length(kxp);
end
SINOamp = single(ones(Nxp,Nxp,nproj));
SINOph = single(zeros(Nxp,Nxp,nproj));
disp('Extracting projection spectra ...')
if nproj<10;  jump = 1;
else          jump = 10; end
for j = 1:nproj
% for j = 90:90 % @@@
    
    if mod(j,jump)==0 || j==nproj; disp([num2str(j) '/' num2str(nproj)]); end
    
	% coords of Ewald spheres in k (illumination) and K=k-k0 (grating) spaces;
	% (Kx,Ky,Kz)=(0,0,0) corresponds to the power peak in Fp at (kx,ky)=(0,0)
	if strcmp(geometry,'fixed')	
		Kx = kxx;
		Ky = kyy;
		if Ewald % kzz==0 outside shifted sphere
            kzz = real(sqrt(kn^2 - (kxx + kxs(j)).^2 - (kyy + kys(j)).^2 )); % in k space
			Kz = kzz - kz(j); % in K=k-k0
            domain = (kzz-4*dkP>0)&(kzz-sqrt(kn^2-(NA/lambda)^2)>0);
		else % infinite-radius Ewald sphere
			kzz = ones(size(kxx)); % just to make kzz>0 (marks domain)
			Kz = kzz - 1;
            domain = (kzz>0);
		end
	elseif strcmp(geometry,'facing')
		Kx = kxx;
		Ky = kyy;
		if Ewald % kzz==0 outside centered sphere
			kzz = real(sqrt(kn^2 - kxx.^2 - kyy.^2)); % in k (projection)
			Kz = kzz - kn; % in K=k-k0
            domain = (kzz-4*dkP>0)&(kzz-sqrt(kn^2-(NA/lambda)^2)>0);
		else % infinite-radius Ewald sphere
			kzz = ones(size(kxx)); % just to make kzz>0 (marks domain)
			Kz = kzz - 1;
            domain = (kzz>0);
		end
		if ~isempty(thetay) % full-angle object rotation
			% rotate the c.s. with left-hand rule (around Ky), so that the data is rotated
			% by theta_x with right-hand rule (around Y):
			Kxr = Kz*sin(thetay(j)) + Kx*cos(thetay(j));
			Kzr = Kz*cos(thetay(j)) - Kx*sin(thetay(j));
			Kx = Kxr;
			Kz = Kzr;
		else % arbitrary angles
			fid = atan2d(kyp(j), kxp(j));
			thd = real(asind( min( 1, sqrt(kxp(j)^2+kyp(j)^2)/(n_imm/lambda) ) ));
			if abs(thd-pi/2)<1e-12; thd=pi/2-1e-12; end % avoid singularity
			RM = SpinCalc('EA123toDCM',-[0 0 fid],1e-6,false) * ...
			     SpinCalc('EA123toDCM',-[0 thd 0],1e-6,false); % theta-POLE
			K_XYZ = RM * [Kx(:).'; Ky(:).'; Kz(:).'];
			Kx = reshape(K_XYZ(1,:), NF, NF);
			Ky = reshape(K_XYZ(2,:), NF, NF);
			Kz = reshape(K_XYZ(3,:), NF, NF);
		end
    end
    
    % discard evanescent waves
    iiF = find(domain);
    % floating point voxel coords
    Kx = Kx/dkP+ceil((NKP+1)/2); % NKP/2+1 corresponds to position of the DC term
    Ky = Ky/dkP+ceil((NKP+1)/2);
    Kz = Kz/dkP+ceil((NKP+1)/2);
    % discard indices of frequencies exceeding the K-cube
    iiF( Kx(iiF)<1 | Kx(iiF)>NKP ) = [];
    iiF( Ky(iiF)<1 | Ky(iiF)>NKP ) = [];
    iiF( Kz(iiF)<1 | Kz(iiF)>NKP ) = [];
    % 1D column vectors for single-subscript addressing
    Kx = Kx(iiF);
    Ky = Ky(iiF);
    Kz = Kz(iiF);
    % nearest neighbor voxels
	iiK = (round(Kz)-1)*NKP^2 + (round(Ky)-1)*NKP + round(Kx);
	if showplots %j==1
	% -------------------------------------------------------------------------
	figure(20) % %(20)!!!!!!!!!!!!!!!!!!!!
	tmp = KO(iiK);
	KO(iiK) = max(KO(:));
	imagesc(log10(lg+abs(squeeze(KO(:,end/2+1,:)).'))); axis image;
% 	set(gcf,'Position',[0 0 1920 1080])
% 	subplot(211); imagesc(log10(lg+abs(squeeze(KO(:,end/2+1,:)).'))); axis image;
% 				  xlabel('Kx'); ylabel('Kz')
% 	subplot(212); imagesc(log10(lg+abs(squeeze(KO(:,:,end/2+1)).'))); axis image;
% 				  xlabel('Kx'); ylabel('Ky')
% 				  colormap jet
	drawnow
	KO(iiK) = tmp;
	% -------------------------------------------------------------------------
	end
	
    % projection spectrum from K space	
    if strcmp(sinterp,'map')
        Fp = single(zeros(NF,NF));
        Fp(iiF) = KO(iiK); % nearest neighbor as in reconstruction  
    else
    %     Fp = single(zeros(NF,NF));
    %     Fp(iiF) = interp3(KO,Ky,Kx,Kz,sinterp); % MATLAB interp3
        Fp_re = single(zeros(NF,NF));
        Fp_im = single(zeros(NF,NF));
        Fp_re(iiF) = interp3(real(KO),Ky,Kx,Kz,sinterp);
        Fp_im(iiF) = interp3(imag(KO),Ky,Kx,Kz,sinterp); % GNU Octave interp3
        Fp=Fp_re+1i*Fp_im;
    end
    if Ewald % Fourier Diffraction Theorem: FFT3{o(x,y,z)}(Ewald) -> FFT2{proj(x,y)}
        if any(kzz<0)
            error('kzz<0 not allowed in FDT formula.')
        end
        kz1 = kzz;
        kz1(~domain) = 1; % do not divide by zero
        Fp = Fp./kz1/(+1i.*(2*pi)*2); % why does + work?
	end
	if showplots
    % ---------------------------------------------------------------------
    figure(21)
	imagesc(log10(lg+abs(Fp.'))); axis image; xlabel('x'); ylabel('y');
% 	imagesc(angle(Fp.')); axis image; xlabel('x'); ylabel('y');% @@@
    title(sprintf('Fp(dkP=%g) %d/%d',dkP,j,nproj))
	colormap jet
	drawnow
    % ---------------------------------------------------------------------
	end
	
% 		Fp = KO_xy; % @@@
    % retrieved projection
    if objshift
        Up = fftshift(ifft2(ifftshift(Fp)))/dx^2;
    else
        Up = ifft2(ifftshift(Fp))/dx^2;
	end
% 	if showplots
%     % ---------------------------------------------------------------------
% 	figure(210)
% 	imagesc(abs(Up'));axis image;xlabel('x');ylabel('y');
% 	title('abs(Up)zpc');colormap jet
%     drawnow
%     % ---------------------------------------------------------------------
% 	end
        
    % reverse zpc (drop signal padding)
    Up = Up((NF-Nx)/2+(1:Nx), (NF-Nx)/2+(1:Nx));
	if objshift
		Fp = fftshift(fft2(ifftshift(Up)))*dx^2;
	else
		Fp = fftshift(fft2(Up))*dx^2;
	end
	if showplots
    % ---------------------------------------------------------------------
%     figure(220);
%     imagesc(abs(Up'));axis image;xlabel('x');ylabel('y');
%     title('abs(Up)');colormap jet
%     drawnow
    figure(22)
	imagesc(log10(lg+abs(Fp.'))); axis image; xlabel('x'); ylabel('y');
    title(sprintf('Fp(dko) %d/%d',j,nproj))
	colormap jet
	drawnow
    % ---------------------------------------------------------------------
    end	
    
    % final size (drop spectrum padding)
	if strcmp(geometry,'fixed') && savemem==0 % center-NA, crop/pad, cirshift-back-NA
		Fp = circshift(Fp, [round(kxs(j)/dk) round(kys(j)/dk)]); % center NA
		Fp = subarray_yx(Fp, [Nxp Nxp], [], 'replicate'); % 'replicate' smudges in gaps due to NKP/zpc~Nxp
		Fp = circshift(Fp, -[round(kxs(j)/dk) round(kys(j)/dk)]); % shift back (wrap)
% 	elseif strcmp(geometry,'facing')
	else % savemem=1: 'fixed'-nothing to wrap anyway, like with 'facing'
		Fp = subarray_yx(Fp, [Nxp Nxp], [], 0);
	end	
    if showplots
    % ---------------------------------------------------------------------
    figure(33)
    imagesc(log10(lg+abs(Fp.'))); axis image; xlabel('x'); ylabel('y');
    title(sprintf('Fp(dko,crop) %d/%d',j,nproj))
	colormap jet
	drawnow	
    % ---------------------------------------------------------------------
	end
	% final projection
	if objshift
        Up = single(fftshift(ifft2(ifftshift(Fp)))/dxp^2);
    else
        Up = single(ifft2(ifftshift(Fp))/dxp^2);
	end
	
	
	if showplots
    % ---------------------------------------------------------------------
    figure(44)
    subplot(221); imagesc(abs(Up.')); title('abs(Up)');  axis image; xlabel('x'); ylabel('y'); colorbar
    subplot(222); imagesc(angle(Up.')); title('angle(Up)'); axis image; xlabel('x'); ylabel('y');
%     subplot(223); plot(abs(Fp(:,end/2+1)), 'b'); axis tight; title('abs(Fp)')
    subplot(223); imagesc(log10(lg+abs(Fp.'))); axis image; title('Fp')
    subplot(224); plot(real(Fp(:,end/2+1)), 'r'); axis tight; title('real(Fp)')
	colormap jet
	drawnow
%	waitforbuttonpress
    % ---------------------------------------------------------------------
	end
	% encoded projection fields --> physical complex amplitudes
    if strcmp(Approx,'PhaseRay')
        SINOph(:,:,j) = abs(Up); % integrated phase
    elseif strcmp(Approx,'Born')
        SINOamp(:,:,j) = abs(Up+1);
        SINOph(:,:,j) = Miguel_2D_unwrapper(single( angle(Up+1) )); % U0=1
    elseif strcmp(Approx,'Rytov')
        amp = abs(exp(Up)); qE = norm(amp(:))/Nxp; % ~ square root of energy quotient (divisor for amp)
        amp = amp/qE;       %qE = norm(amp(:))/Nxp; % correct projection energy flux to unity (& check)
        ph = Miguel_2D_unwrapper(single( angle(exp(Up)) ));
		ph_sign = sign(sum(ph(:)));
		if ph_sign~=obj_sign
			ph = ph + obj_sign*2*pi; % sometimes Miguel assigns ph=0 to object center
        end        
        SINOamp(:,:,j) = amp;
		SINOph(:,:,j) = ph;
    else
        error('Wrong option `Approx`');
	end
    	
	if showplots && figs_autoarrange && j==1
		autoArrangeFigures()
	end
	
end


