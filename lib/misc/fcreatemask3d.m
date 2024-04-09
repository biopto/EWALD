function[Y] = fcreatemask3d(X)

% X         - input 3D array

X(X<0)=0;
X = mat2gray(X);

level = graythresh(X(:));
X(X>level) = 1; X(X<=level) = 0;
mask = +X;

% dilate mask
se = strel('sphere',3);
mask = imdilate(mask,se);

n = 3;
H = 1/(n^3)*ones(n,n,n);
Y = single(convn(+mask,H,'same'));

end