% Resize n-dimensional array to nsize==[N1 N2 ... Nn]
% Works best when all dimensions are rescaled at the same time, otherwise one gets the 'pixelosis effect'
% (see sino3DCropRes.m for sinogram resampling)
%
% Copyright: 2015, Piotr L. Makowski, p.makowski@mchtr.pw.edu.pl
%

function mat = ndresize(mat, nsize, interp_mode)
% mat           - N-dimansional numerical array
% nsize         - new size list
% interp_mode   - {'nearest', 'linear', 'cubic', 'spline'}

if length(nsize)~=ndims(mat)
    error('`nsize` should be an array of ndims(mat) elements')
end
osize = size(mat);
nsize = round(nsize);

if sum(size(mat) ~= nsize)    

    n = ndims(mat);
    qnc = cell(1,n);
    for j=1:n
        qnc{j} = single( (0:nsize(j)-1)/max(1,(nsize(j)-1)) * (osize(j)-1) + 1 );
    end
    [qqn{1:n}] = ndgrid(qnc{:});

	mat_re = interpn(real(mat), qqn{:}, interp_mode);
    mat = 1i*interpn(imag(mat), qqn{:}, interp_mode);
	mat = mat + mat_re;% Octave
    
end

