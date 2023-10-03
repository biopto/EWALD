disp('Saving reconstruction...')
    
    FileName = 'REC[';
    if Masking
        FileName = strcat(FileName,'GPSC]');
    elseif nGPi == 0
        FileName = strcat(FileName,'DI]');
    else
        FileName = strcat(FileName,'GP]');
    end
    
    % temporarly added due to error for not empty n_obj 
    % (temporarly forced to be empty within defaults.m script) 
    % and Kspace_padding > NaN:
    % if isempty(n_obj); sinogram = strrep(sinogram, '_nobj', ''); end
    FileName = strcat(FileName,'-' ,sinogram(1:end-4));
    %FileName = sprintf('%s-%s_Kpad%.2f', FileName, sinogram(1:end-4), Kspace_padding);
%     if exist('rec_mode','var') && contains(FileName,'OCDT')
%         if strcmp(rec_mode,'OCT') || strcmp(rec_mode,'ODT')
%             FileName = strrep(FileName,'OCDT',rec_mode);
%         end
%     end
    if sinogram_subset_factor > 1
         FileName = sprintf('%s_(subs%d)', FileName, length(sino_params(1,:))); 
    end
    % FileName = sprintf('%s_Kpad%.2f', FileName, Kspace_padding); 
    
    str_file = [PathName FileName '.mat'];
    disp(str_file)
    
    solverParams = {'time',sim_time;'projection_crop_factor',projection_crop_factor;...
                    'projection_downsample_factor',projection_downsample_factor;...
                    'N_SINO_y',N_SINO_y;'N_projection_padded',N_projection_padded;...
                    'projection_padding_xy',projection_padding_xy;...
                    'Kspace_oversampling_z',Kspace_oversampling_z;'Kspace_padding',Kspace_padding;...
                    'nGPi',nGPi;'relaxGP',relaxGP;'relaxM',relaxM;...
                    'do_NNC',do_NNC;...
                    'RMAEtab',RMAEtab;'RRMSEtab',RRMSEtab;...
                    'RMADtab',RMADtab;'RRMSDtab',RRMSDtab;'epsi',epsi};
      % save(str_file, 'RECON', 'dxo', 'n_imm', 'oct_ind', 'solverParams')
      save(str_file, 'RECON', 'dx', 'n_imm', 'sino_params', 'solverParams','-v7.3')