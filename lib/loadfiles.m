[sinogram,PathName] = uigetfile('*.mat','Select sinogram');
if sinogram~=0
    load([PathName sinogram]); % (x,y,proj)  
else; error('You have to load a sinogram first.');
end    

% Modify names of parameters:
if ~exist('NA','var') && exist('NA_','var')
    NA = NA_; clear NA_
end
if ~exist('geometry','var') && exist('geometry_','var')
    geometry = geometry_; clear geometry_
end

% Change the name of immersion refractive index parameter:
if exist('n_immersion','var') 
    n_imm = n_immersion; clear n_immersion;
else; error('Immersion medium refractive index value not found.'); 
end
