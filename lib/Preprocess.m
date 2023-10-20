function [SINOamp_reduced,SINOph_reduced, sino_params, dx, ...
    projection_padding_xy, N_projection_padded, ...
    x_tv, n_obj_res, plots, Fpmask, projection_downsample_factor, N_SINO_y] = ...
			Preprocess(SINOamp, SINOph, geometry, sino_params, resample_projections, dx, projection_crop_factor, Fpmask,...
            n_obj, thetay, sinogram_subset_factor, n_imm, Masking, N_CP, nCPi,...
            projection_padding_xy, plots, do_NNC) 
          
%% Data preprocessing: cropping and resampling input data
if ~exist('SINOamp','var')
    SINOamp = single(ones(size(SINOph))); % (x,y,proj)
end

%% reduce SINO (select subset of projections)

if isempty(thetay)
	[sino_params, ind] = ...
                        reduce_sources_2pi(sino_params, sinogram_subset_factor);
    rx_reduced = sino_params(1,:); ry_reduced = sino_params(2,:);
    sino_params(2,:) = n_imm./sino_params(3,:) .* sino_params(2,:);
	sino_params(1,:) = n_imm./sino_params(3,:) .* sino_params(1,:);
else
	ind = 1 : sinogram_subset_factor : length(thetay);
	thetay = thetay(ind);
end
SINOph_reduced = SINOph(:,:,ind);
SINOamp_reduced = SINOamp(:,:,ind);
clear SINOph SINOamp
if ~isempty(Fpmask)
    Fpmask = Fpmask(:,:,ind);
end

%% Resample projections
switch geometry
	case 'fixed'
        dx_max = min( sino_params(3,:)./(4*sino_params(4,:)) );
	case 'facing'
        dx_max = min( sino_params(3,:)./(2*sino_params(4,:)) );
        % always optimal, but downsampling to this level
        % for 'fixed' must involve spectral circshifting
        % to avoid data loss; this should be done in sinogram
        % preparation, not worth the effort at this point
end

if(resample_projections && dx_max>dx) 
    projection_downsample_factor = dx_max/dx; 
else
    projection_downsample_factor = 1; 
end

if(projection_downsample_factor>1 || projection_crop_factor<1)
	SINOamp_reduced = sino3DCropRes(SINOamp_reduced, projection_crop_factor, projection_downsample_factor);
	SINOph_reduced = sino3DCropRes(SINOph_reduced, projection_crop_factor, projection_downsample_factor);
    if ~isempty(Fpmask)
        Fpmask = sino3DCropRes(Fpmask, projection_crop_factor, projection_downsample_factor);
    end
	dx = dx*projection_downsample_factor;
end

%% match size and resolution of phantom object with projections
N_SINO_y = size(SINOph_reduced,1);
if ~isempty(n_obj)% (x,y,z)
	n_obj_res = ndresize(n_obj, [N_SINO_y N_SINO_y N_SINO_y], 'linear'); % match resolution with original SINO
	n_obj_res = ndcrop(n_obj_res, size(n_obj_res)*projection_crop_factor); % croppng before resampling (sino3DCropRes)
	n_obj_res = ndresize(n_obj_res, [N_SINO_y N_SINO_y N_SINO_y], 'linear'); % match with SINO
else
	n_obj_res = [];
end


%% Stage 1 - Generate object support
if Masking
    % Generate illumination versors from illumination directions for 1st stage of GTVIC
    if isempty(thetay)
        vectors = rays_to_vectors_new(rx_reduced,ry_reduced, geometry);
    end
    
	% ASTRA sinogram convention for Chambolle-Pock: (x,proj,y)
	SINOamp_reduced = permute(SINOamp_reduced,[1 3 2]);
	SINOph_reduced = permute(SINOph_reduced,[1 3 2]);
		
    % ASTRA parameters
    vol_geom = astra_create_vol_geom(N_CP, N_CP, N_CP);
    proj_geom = astra_create_proj_geom('parallel3d_vec', N_CP, N_CP, vectors);

    % RESIZING SINOGRAM FOR QUICK CALCULATION OF CP
    CP_downsampling = N_SINO_y/N_CP;

    SINOph_CP = zeros(N_CP, size(SINOph_reduced,2), N_CP); % (x,proj,y)
    ny=N_CP; nx=N_CP; %% desired output dimensions
    [x, y]=ndgrid(linspace(1,size(SINOph_reduced,1),nx),...
                  linspace(1,size(SINOph_reduced,3),ny));
    for ii=1:size(SINOph_reduced,2) 
        SINOph_CP(:,ii,:) = ...
            interp2(...
                reshape(SINOph_reduced(:,ii,:), [size(SINOph_reduced,1) size(SINOph_reduced,3)]) ...
            ,y,x,'cubic' );
    end
    % phase normalization for SIRT
    SINOph_CP = SINOph_CP./(dx*CP_downsampling)./(2*pi./sino_params(3,:));

    %% Stage 1.1 - Generation of an object suppport	
    A = opTomo('cuda', proj_geom, vol_geom);
    TV3D = opTV3D(N_CP) ;
    starting_point = zeros(N_CP,N_CP,N_CP);
    scale = 0.05; %0.05 for demon
    x_tv = chambolle_pock(...
                    scale*A, TV3D, (1-2*(do_NNC<0))*scale*SINOph_CP(:), nCPi, 0.002, true, starting_point(:));
    clear SINOph_CP starting_point TV3D A
    x_tv = (1-2*(do_NNC<0))*x_tv;
    x_tv = single(reshape(x_tv,[N_CP N_CP N_CP]));

	
    %% Resample object support to match projection size
	x_tv = ndresize(x_tv, [N_SINO_y N_SINO_y N_SINO_y], 'linear');
	
	%% Back to FDT sinogram convention: (x,proj,y)->(x,y,proj)
	SINOamp_reduced = permute(SINOamp_reduced,[1 3 2]);
	SINOph_reduced = permute(SINOph_reduced,[1 3 2]);

else
    x_tv = [];
end 

%% Stage 2 - ODT Reconstruction
disp('[Reconstruction]')
N_projection_padded = round(projection_padding_xy*N_SINO_y/2)*2;% padded projection size

plots = upper(num2str(plots));

end