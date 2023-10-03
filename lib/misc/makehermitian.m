function [herm] = makehermitian(mat)

% if ndims(mat)~=3
%     error('Processed matrix must have 3 dimensions!')
% elseif ~isequal(size(mat,1),size(mat,2),size(mat,3))
%     error('Processed matrix must have equal dimensions!')
% end

n = size(mat,1);

%2D
% mat_conj = zeros(n,n);
% mat_conj(n:-1:1,n:-1:1) = mat(1:n,1:n);
% mat_conj = circshift(mat_conj,1,1);
% mat_conj = circshift(mat_conj,1,2);
%3D
mat_conj = zeros(n,n,n);
mat_conj(n:-1:1,n:-1:1,n:-1:1) = mat(1:n,1:n,1:n);
mat_conj = circshift(mat_conj,1,1);
mat_conj = circshift(mat_conj,1,2);
mat_conj = circshift(mat_conj,1,3);

mat_conj = conj(mat_conj);

herm = mat_conj + mat;

ind = find(abs(mat_conj)>0 & abs(mat)>0);
herm(ind) = herm(ind)/2;
