% 1. Load old-type ocdt sinogram
% 2. Run the following lines:
oct_ind = find(rayXY(1,:) > max_range/2);
rayXY(:,oct_ind) = rayXY(:,oct_ind) - offset;

sino_params(1:2,:) = [rayXY]; % directions
sino_params(3,  :) = lambda_all; % wavelengths
sino_params(4,  :) = ones(1,size(rayXY,2)).*NA;
sino_params(5,  :) = ones(1,size(rayXY,2)).*2;

clear rayXY octNA oct_lambda NA lambda lambda_all lambda_vec offset oct_ind Nx nproj max_range Ewald single_c

% 3. Save the resulting data as the new-type sinogram 