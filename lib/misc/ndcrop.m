% Crop n-dimensional array to nsize==[N1 N2 ... Nn]
%
% Copyright: 2015, Piotr L. Makowski, p.makowski@mchtr.pw.edu.pl
%

function [mat] = ndcrop(mat, nsize)
% mat		- N-dimensional array
% nsize		= [n1,n2,...,nN] new size


if ndims(mat)~=length(nsize)
    error('`nsize` should be an array of ndims(mat) elements')
end
osize = size(mat);
nsize = round(nsize);

% disp('N-D cropping ...')
% disp(num2str(osize))
% disp(num2str(nsize))

n = ndims(mat);
args = cell(1,n);
for j=1:n
    if (nsize(j)<1 || nsize(j)>osize(j))
		error('nsize=%g must be <= osize=%g',nsize(j),osize(j));
	end
    args{j} = (1:nsize(j)) + round( (osize(j)-nsize(j))/2 );
end

mat = mat(args{:});
