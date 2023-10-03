% Reduction of number of sources/projections
% author: Wojciech Krauze and Pawel Ossowski, WUT
% 
%%%%%%%%%%%% INPUT
% sino_params  -  sinogram parameters including kx,ky (corrdinate system center is 0,0 )
% proj_step    -  which projections should be taken
%                 e.g. 2 means that every 2nd projection will be taken
%%%%%%%%%%%% OUTPUT
% sino_params  - updated sinogram parameters
% ind          - vector which specifies how a sinogram should be
%                rearranged to be compatible with updated parameters

function [sino_params, ind] = ...
                          reduce_sources_2pi(sino_params, proj_step)

kx = reshape(sino_params(1,:), [], 1);
ky = reshape(sino_params(2,:), [], 1);

fi = atan2(-ky,kx); % -ky to align angular order with ASTRA order (center-left -> down)

[~,ind] = sort(fi) ; % fi_sort = (fi(ind))
% ind = circshift(ind, [1 0]); % start from fi=2pi (closest to X axis)
sino_params = sino_params(:,ind);

% lambda_all_sort = lambda_all(ind);

% lambda_vec = unique(lambda_all);
% ind_centr = find(lambda_vec==median(lambda_vec)); % central wavelength index in lambda_vec vector

% vec_proj_odt = vec_proj(sino_params(5,:)==1);
% vec_proj_oct = vec_proj(sino_params(5,:)==2);

% if ~isnan(wavelengths_odt) && (wavelengths_odt>=0 && wavelengths_odt<=length(lambda_vec))
%     % Reduce the number of wavelengths due to `wavelengths_odt` parameter (around central wavelength value):
%     lambda_vec_odt = lambda_vec( ind_centr-floor(wavelengths_odt/2) : ind_centr+(wavelengths_odt-1)-floor(wavelengths_odt/2) );
%     vec_proj_odt = vec_proj_odt(ismember(lambda_all_sort(odt_ind), lambda_vec_odt)); 
% end
% if ~isnan(wavelengths_oct) && (wavelengths_oct>=0 && wavelengths_oct<=length(lambda_vec))
%     lambda_vec_oct = lambda_vec( ind_centr-floor(wavelengths_oct/2) : ind_centr+(wavelengths_oct-1)-floor(wavelengths_oct/2) );   
%     vec_proj_oct = vec_proj_oct(ismember(lambda_all_sort(oct_ind), lambda_vec_oct));
% end
% vec_proj = sort([vec_proj_odt vec_proj_oct]);

% error('Parameter value out of range of [0 %d]. Change `wavelengths_odt` parameter value.', length(lambda_vec));
% vec_proj = vec_proj(ismember(lambda_all, lambda_vec_red));     

% projection number reduction (if proj_step>1):
selected_proj = 1:proj_step:length(kx);
sino_params = sino_params(:,selected_proj);
ind = ind(selected_proj);

end