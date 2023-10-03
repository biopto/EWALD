% Saves reconstruction output to .png and reoconstruction data to .mat.
% Run this manually after axecuting 'clear' and loading previously saved
% Rec[...].mat to recreate plots, change names, etc.
%__________________________________________________________________________
% 03.12.2018
% - automatically flip dnlims_man for do_NNC<0
% - when dnlims_man=[] choose dnmax according to OBJ, if available
% - nlims~2*clims consistently
% 01.10.2018
% - fix bugs with center_cube_3d
% 31.07.2018
% - remove: rng, Nxo
% - fix bugs with ndlims_man and clims
% - add automatic centering of object in REC inside a cube of size Nc
% - add projection_padding_xy to saved parameters (spectrum oversampling factor)
% - remove qcrop, REC cropping controlled by projection_padding_xy
% 28.02.2018
% - correct str2 bug for MASKING=-1
% - improve processing of n_obj for more accurate comparison with rec_complex
% 05.12.2017
% - fix plotting limits bug for negative REC-n0
% 14.11.2017
% - detect manual run by absence of SINOph_reduced, no need to clear between runs
% - save only slices of KO and KOi instead of entire clubes
% - simplify automatic nlims, fix bug for nlims_man(1)<0
% 25.10.2017
% - fix OBJ size bug in case of projection padding: N_projection_padded=/=N_SINO_y
% - bring qcrop option here (REC crop coefficient for display only)
% - add to saved params: projection_crop_factor, projection_downsample_factor, N_SINO_y
% 13.09.2017
% - fix OBJ size for projection_crop_factor<1 (cropped projections)
% - move all parameters to beginning
% - display OBJ cross-sections as area plots
% 10.07.2017
% - add options 'saverecs', 'makeplots'
% - add 'qcrop' to saved params (manual crop defined in RecGTVIC)
% 30.06.2017
% - do not save PathName; when run manually, use current path (pwd)
% - add N_SINO_y,N_projection_padded to variables saved in 'Rec' file to account for N_projection_padded in manual run
% - add section to display OBJ instead of REC (commented out by default)

% comment = sprintf('_(%s-x%gz%g)RMAD_relaxLaR',Interp,N_projection_padded/N_SINO_y,Kspace_oversampling_z);
% comment = '';
% if(Masking && ~GTVIC); MASKING = -Masking; else; MASKING = Masking; end
if Masking; MASKING = -Masking; else; MASKING = Masking; end 
N_Kspace_xy_padded = size(RECON,1);
N_Kspace_z_padded_upsampled = size(RECON,3);
Nc = NaN;% cropping cube
%%%%%%%%%%%%%%%%%%%%%%%

is_n0 = true;%true % whether n_obj already holds absolute refractive index
show_n0 = false; % what to display

tick = 0.02; % for automatic RI range
dnlims_man = [];%[.0 .05];%[1.48 1.52]; % manual clims for REC slices (set [] for clims='auto')
if do_NNC<0 && ~isempty(dnlims_man); dnlims_man(1) = -dnlims_man(2); dnlims_man(2) = .0; end

% logadd = 0.0001; % 0.1 % offset for logarithmic plots
logadd = 0.1;
wp = 1; % plot line width
opaq = .2; % opaqueness of area plots
revcolors = false; % reverse jet colormap (for rec_complex < n_imm)

% cross-section coords in [N_Kspace_xy_padded N_Kspace_xy_padded N_Kspace_xy_padded] cube
cx0 = round(N_Kspace_xy_padded/projection_padding_xy/2)+1;
cy0 = round(N_Kspace_xy_padded/projection_padding_xy/2)+1;
cz0 = round(N_Kspace_xy_padded/projection_padding_xy/2)+0;
% % manual
% cx0 = 202;%202(Meth)
% cy0 = 194;%194(Meth)
% cz0 = 150;%150(Meth)
% % auto
% [cx0,cy0,cz0] = center_mass_3d(max(0, do_NNC*real(ndcrop(rec_complex,round([N_Kspace_xy_padded N_Kspace_xy_padded N_Kspace_xy_padded]/projection_padding_xy/2)*2)-is_n0*n_imm)));
% % take a smaller cube centered aroud (cx,cy,cz)
% Nc = 200;%100
%%%%%%%%%%%%%%%%%%%%%%%


if all(sino_params(5,:)==2) 
    Kz_slice = 50;
else; Kz_slice = size(KO,3)/2+1; 
end


% Collect reconstruction info (x,y,z)
N_projection_padded = round(N_projection_padded/2)*2;
if exist('SINOph_reduced','var') % when called by RecGTVIC
	KO_xz = squeeze(KO(:,end/2+1,:));
	KOi_xz = squeeze(KOi(:,end/2+1,:));
% 	KO_xy = squeeze(KO(:,:,end/2+1));
% 	KOi_xy = squeeze(KOi(:,:,end/2+1));
    KO_xy = squeeze(KO(:,:,Kz_slice));
	KOi_xy = squeeze(KOi(:,:,Kz_slice));
    
% 	if ~exist('SINOamp_reduced','var') || phase_only
% 		strPH = '-PH';
% 	else
% 		strPH = '';
% 	end
	%str1 = [sinogram(1:end-4) strPH sprintf('(..%g)',sinogram_subset_factor)];
%  	str2 = sprintf('c%.1f-%d.%d-%sE%d%s(i%d)(j%d-i%gi%g-%g)nnc%smsk%s', ...
%  				projection_crop_factor,N_projection_padded,N_Kspace_xy_padded, geometry,Ewald,Approx, ...
%  				nGPi,nPEi,(nPEi>0)*nEi1,(nPEi>0)*nEi2,(nPEi>0)*epc, num2str(do_NNC), num2str(MASKING));
% 	str3 = sprintf('lam%.1fNA%.2f(n%.4f)dx%.3f',lambda*1e3,NA,n_imm,dx_max);
else % when 'Rec' file loaded manually
    PathName = pwd; PathName = [PathName '/']; % current folder
%     str1 = '';
%     str2 = '';
%     str3 = '';
% 	comment='';
% 	if strfind(str2, 'msk1'); MASKING = 1;
% 	elseif strfind(str2, 'msk-1'); MASKING = -1;
% 	else MASKING = 0;end
end

% str1;
% str2 = [str2 comment];
% str3;

% Plot reconstruction data (y,x,z)

% OBJECT PHANTOM (y,x,z)
n0 = n_imm;
if ~isempty(n_obj)
    n_obj_res = ndresize(n_obj, [N_SINO_y N_SINO_y N_SINO_y], 'linear'); % match resolution with original SINO
    n_obj_res = ndcrop(n_obj_res, size(n_obj_res)*projection_crop_factor); % croppng before resampling (sino3DCropRes)
    n_obj_res = ndresize(n_obj_res, size(n_obj_res)./projection_downsample_factor, 'linear'); % resampling (sino3DCropRes)
    n_obj_res = padarray(n_obj_res, [.5 .5 .5].*(N_projection_padded-N_SINO_y), is_n0*n0); % projection padding
    OBJ = ndresize(permute(real(n_obj_res), [2 1 3]), [N_Kspace_xy_padded N_Kspace_xy_padded N_Kspace_xy_padded], 'linear') ...
            -(~is_n0)*n0; % (y,x,z)
    OBJ = ndcrop(OBJ, round(size(OBJ)/projection_padding_xy /2)*2);
else
    OBJ = [];
end
% OBJ,REC cropping
[RECON, cy,cx,cz] = center_cube_3d(RECON, cy0,cx0,cz0, Nc, is_n0*n0);
[OBJ, ~,~,~] = center_cube_3d(OBJ, cy0,cx0,cz0, Nc, is_n0*n0);

if ~isempty(n_obj)		
    tmp_rec = OBJ(...
                round(0.3*size(RECON,1)) : round(0.7*size(RECON,1)),...
                round(0.3*size(RECON,1)) : round(0.7*size(RECON,1)),...
                round(0.3*size(RECON,1)) : round(0.7*size(RECON,1)))-is_n0*n0;
else
    tmp_rec = RECON(...
                round(0.3*size(RECON,1)) : round(0.7*size(RECON,1)),...
                round(0.3*size(RECON,1)) : round(0.7*size(RECON,1)),...
                round(0.3*size(RECON,1)) : round(0.7*size(RECON,1)))-is_n0*n0;
end
nm = max(abs(tmp_rec(:)));
clear tmp_rec
dnmax = max([tick, tick*ceil((nm-tick/2)/tick)]); % estimate of max(delta_n)
if isempty(dnlims_man)
    nlims = [-dnmax dnmax]*1.5 + show_n0*n0;
    clims = 'auto';
else
    nlims = [-max(abs(dnlims_man-show_n0*n0)) max(abs(dnlims_man-show_n0*n0))]*2.0 + show_n0*n0;
    clims = dnlims_man;
end

%  slices
corr_n0 = +(~is_n0&show_n0)*n0 -(~show_n0&is_n0)*n0;% correction for n0
rec_zx = squeeze(RECON(cy,:,:) + corr_n0).';%(prim: make X horizontal)
rec_yx = squeeze(RECON(:,:,cz) + corr_n0);
if ~isempty(OBJ) % calculate Quality Index for one slice
    obj_zx = squeeze(OBJ(cy,:,:) + corr_n0).'; %  one slice (prim: make X horizontal)
    diff_zx = rec_zx - obj_zx;
    RMS_zx = rms(diff_zx(:));
    QI_zx = Q_QualityIndex(obj_zx,rec_zx);

    obj_yx = squeeze(OBJ(:,:,cz) + corr_n0);
    diff_yx = rec_yx - obj_yx;
    RMS_yx = rms(diff_yx(:));
    QI_yx = Q_QualityIndex(obj_yx,rec_yx);

    wp = 2.5;
end


if Masking; str4 = 'GPSC'; elseif nGPi==0; str4 = 'DI';
else; str4 = 'GP'; end

for ii=1:2 % repeat to correct subplot sizes
    figure(100)
    clf

    % 1st chart row
    subplot(241); imagesc(log10(logadd +abs(KO_xz'))); axis square; %axis image
                  xlabel('Kx'); ylabel('Kz')
                  %title([strrep(str1,'_','-'), '                                   '])
                  title 'KO\_xz (initial KO)'
    subplot(242); imagesc(log10(logadd +abs(KOi_xz'))); axis square; %axis image
                  xlabel('Kx'); ylabel('Kz')
                  %title(strrep(str2,'_','-'))
                  title(sprintf('KOi\\_xz - final rec (%s)', str4));
    [nz, nx] = size(rec_zx);
    subplot(243); imagesc(real(rec_zx)); axis image; %colorbar
                  hold on
                  plot([1,nx],[cz,cz],':','Color',[.5 .5 .5],'LineWidth',1)
                  hold off
                  xlabel('x'); ylabel('z')
                  %title(sprintf('                    cx0=%d,  cy0=%d,  cz0=%d',cx0,cy0,cz0))
                  title 'RECON x-z'
                  caxis(clims)
    subplot(244); hold on
                  if ~isempty(OBJ)
                  plot(real(obj_zx(cz,:)), 'b','LineWidth',1)
                  plot(real(obj_zx(:,cx)), 'r','LineWidth',1)
    % 				  h1=area(real(obj_zx(cz,:)),'ShowBaseLine','off','LineStyle','none','EdgeColor',[0 0 1],'LineWidth',1);
    % 				  h2=area(real(obj_zx(:,cx)),'ShowBaseLine','off', 'LineStyle','none','EdgeColor',[1 0 0],'LineWidth',1);
    % 				  set(h1,'FaceAlpha',opaq,'FaceColor',[0 0 1]);
    % 				  set(h2,'FaceAlpha',opaq,'FaceColor',[1 0 0]);
                  end
                  fig1=plot(real(rec_zx(cz,:)), 'b','LineWidth',wp);
                  fig2=plot(real(rec_zx(:,cx)), 'r','LineWidth',wp);
    % 				  fig1=area(real(rec_zx(cz,:)), 'LineStyle','none');
    % 				  fig2=area(real(rec_zx(:,cx)), 'LineStyle','none');
    % 				  set(fig1,'FaceAlpha',opaq,'FaceColor',[0 0 1]);
    % 				  set(fig2,'FaceAlpha',opaq,'FaceColor',[1 0 0]);
                  hold off
                  box on; grid on; axis tight; legend([fig1 fig2], {'x','z'}, 'Location', 'SouthEast'); 
    %               box on; grid on; axis tight; legend([fig1 fig2], {'x','z'}, 'Location', 'SouthEast');
                  title 'REC x,z profiles'
                  ylim(nlims);
    % 			  ylim([-1 1]*0.07); 
    if ~isempty(OBJ)
                  title(sprintf('REC x,z profiles; RMS(xz)=%.2fe-3, QI(xz)=%.2f%%',RMS_zx*1e3,QI_zx*100))
    end

    % 2nd chart row
    subplot(245); imagesc(log10(logadd +abs(KO_xy'))); axis image
                  xlabel('Kx'); ylabel('Ky'); colormap jet
                  title 'KO\_xy (initial KO)'
    subplot(246); imagesc(log10(logadd +abs(KOi_xy'))); axis image
                  xlabel('Kx'); ylabel('Ky'); colormap jet
                  %title(strrep(str3,'_','-'))
                  title(sprintf('KOi\\_xy - final rec (%s)', str4));
    [ny, nx] = size(rec_yx);
    subplot(247); imagesc(real(rec_yx)); axis image; %colorbar
                  hold on
                  plot([1,nx], [cy,cy],':','Color',[.5 .5 .5],'LineWidth',1)
                  hold off
                  xlabel('x'); ylabel('y')
                  title 'RECON x-y'
    % 				  title(sprintf('cz=%d',cz)) 
                  caxis(clims)
    if ~isempty(OBJ)
                  title(sprintf('REC x-y; MAE=%.3fe-3, MSE=%.3fe-3',RMAEtab(nGPi+1)*1e3,RRMSEtab(nGPi+1)*1e3))
    				  RMAE_REC = norm(RECON(:)-OBJ(:),1)/norm(OBJ(:),1);
    				  RRMSE_REC = norm(RECON(:)-OBJ(:),2)/norm(OBJ(:),2);
    				  title(sprintf('MAE=%.3fe-3, MSE=%.3fe-3',RMAE_REC*1e3,RRMSE_REC*1e3))
    end
    subplot(248); hold on
                  if ~isempty(OBJ)
                  plot(real(obj_yx(cy,:)), 'b', 'LineWidth',1)
                  plot(real(obj_yx(:,cx)), 'r', 'LineWidth',1)
    % 				  h1=area(real(obj_yx(cy,:)),'ShowBaseLine','off','LineStyle','none','EdgeColor',[0 0 1],'LineWidth',1);
    % 				  h2=area(real(obj_yx(:,cx)),'ShowBaseLine','off','LineStyle','none','EdgeColor',[1 0 0],'LineWidth',1);
    % 				  set(h1,'FaceAlpha',opaq,'FaceColor',[0 0 1]);
    % 				  set(h2,'FaceAlpha',opaq,'FaceColor',[1 0 0]);
                  end
                  fig1=plot(real(rec_yx(cy,:)), 'b', 'LineWidth',wp);
                  fig2=plot(real(rec_yx(:,cx)), 'r', 'LineWidth',wp);
    % 				  fig1=area(real(rec_yx(cy,:)), 'LineStyle','none');
    % 				  fig2=area(real(rec_yx(:,cx)), 'LineStyle','none');
    % 				  set(fig1,'FaceAlpha',opaq,'FaceColor',[0 0 1]);
    % 				  set(fig2,'FaceAlpha',opaq,'FaceColor',[1 0 0]);
                  hold off
                  box on; grid on; axis tight; legend([fig1 fig2], {'x','y'}, 'Location', 'SouthEast');
    %               box on; grid on; axis tight; legend([fig1 fig2], {'x','y'}, 'Location', 'SouthEast');
                  title 'REC x,y profiles'
                  ylim(nlims);
    % 			  ylim([-1 1]*0.07);
    if ~isempty(OBJ)
                  title(sprintf('REC x,y profiles; RMS(xy)=%.2fe-3,  QI(xy)=%.2f%%',RMS_yx*1e3,QI_yx*100))
    end

    % progress data
    if nGPi>0
        subplot(245);
                    plot(0:length(RMADtab)-1, RMADtab, 'kx-');hold on
                    plot(0:length(RRMSDtab)-1, RRMSDtab, '+-');hold off
                    grid on
                    ylim([0 5e-5]);xlim([0 length(RMADtab)])
                    legend({'RMAD','RRMSD'},'Location','northeast')
        %             title(sprintf('epsi=(%g, %g, %g),  time=%.1fmin', ...
        %                     (nPEi==0)*epsi+(nPEi>0)*nGPi,(nPEi>0)*nEi1,(nPEi>0)*nEi2,time/60))
    %                     title(sprintf('KO\\_xy; epsi=(%g, %g, %g), %s', ...
    %                             (nPEi==0)*epsi+(nPEi>0)*nGPi,(nPEi>0)*nEi1,(nPEi>0)*nEi2,sim_time))
%                     title(sprintf('epsi=(%g, %g, %g)', ...
%                                  (nPEi==0)*epsi+(nPEi>0)*nGPi,(nPEi>0)*nEi1,(nPEi>0)*nEi2))
        if ~isempty(OBJ)
            subplot(246);
                        plot(0:length(RMAEtab)-1,RMAEtab,'k-o');hold on
                        plot(0:length(RRMSEtab)-1,RRMSEtab,'-^');hold off
                        grid on
                        ylim([0 5e-3]);xlim([0 length(RMAEtab)])
                        legend({'RMAE','RRMSE'},'Location','northeast')
                        title(sprintf('relaxGP=%.2f',str4,relaxGP))
                        % title(sprintf('KO Kx-Ky (%s)', str4));
        end
    end

    pos = [0 0 1920 1080]; set(gcf, 'Position', pos)
end % repeat for better shape
clear str4

if revcolors; colormap(flipud(colormap)); end

% file_png = [PathName 'Plots-REC[' str4 ']-' sinogram(1:end-4) '.png'];
% file_png = sprintf('%sPlots-REC[%s]-%s.png', PathName, str4, sinogram(1:end-4));
% print_plot(gcf, [PathName 'Plots-REC[' str4 ']-' sinogram(1:end-4) '.png'], 100, pos);

if exist('FileName','var') && exist('PathName','var')  
    print_plot(gcf, [PathName 'Plots-' FileName '.png'], 100, pos);
end

% if MASKING
%     print_plot(gcf, [PathName 'P1--' str1 '-' str2 '.png'], 100, pos)
% else
%     print_plot(gcf, [PathName 'p1--' str1 '-' str2 '.png'], 100, pos); 
% end
%  	saveas(gcf,['P1--' str1 '-' str2 '.fig']);

%	figure;imagesc(real(obj_yx));axis image;caxis(clims);title('obj');colormap jet
%	figure;imagesc(real(rec_yx));axis image;caxis(clims);title('rec');colormap jet

