function[Y] = fcreatemask3d(X,type)

% X         - input 3D array
% type      - object type ('adaptive' or 'otsu')

    % S = size(X) ;
    % [x, y, z] = meshgrid(1:S(1),1:S(2),1:S(3)) ;
    % mask_sphere = sqrt( (x-round(S(1)/2)).^2 + (y-round(S(2)/2)).^2 +(z-round(S(3)/2)).^2)<=round(min(S)/2);
    X(X<0)=0;
    X = mat2gray(X);
    if strcmp(type,'adaptive')
        X = adaptthresh(X,'neigh',[15 15 3],'Fore','bright');
        X(X>0.01) = 1; X(X<=0.01) = 0;
        % mask = +X.*mask_sphere;
        mask = +X;
    elseif strcmp(type,'otsu')
        level = graythresh(X(:));
        X(X>level) = 1; X(X<=level) = 0;
        % mask = +X.*mask_sphere;
        mask = +X;
%         dilate mask
        se = strel('sphere',3);
        mask = imdilate(mask,se);
%         mask(:,:,round(S(3)/2-5):round(S(3)/2+5)) = 1;
    else 
        error('Wrong segmentation method (must be ''adaptive'' or ''otsu'')');
    end
    
    n = 3;
    H = 1/(n^3)*ones(n,n,n);
    Y = single(convn(+mask,H,'same'));

end