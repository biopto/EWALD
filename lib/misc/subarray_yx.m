% Copyright: 2018, Piotr L. Makowski, p.makowski@mchtr.pw.edu.pl

function [mat] = subarray_yx(mat, n_yx, shift_yx, padval)
% [mat] = subarray_yx(mat, n_yx, shift_yx, padval)
%
% Extract a shifted sub-array, fill unknown regions with padval
% Convention:
% 	even array/subarray: ...,-2,-1,0,1,...
% 	odd array/subarray:  ...,-1,0,1,...
%
% %%%%%%%%% INPUT
%	mat       - 2D array
%	n_yx      = [noy,nox] output frame size;
%	            set NaN for one of sizes to keep it fixed;
%	            set [] to keep original size intact
% %%%%%%%%% <optional>
%	shift_yx  = [Dy,Dx] frame position shift with respect to
%	            input array center, in pixels
%	padval    - (default=0) padding value/method for padarray()
%
% Copyright: 2016-2018, Piotr L. Makowski, p.makowski@mchtr.pw.edu.pl
%
% -----------------------------------------
% mat = [[1 2 3 4 5 6];  [1 2 3 4 5 6].*10]
% n_yx = [4 4];
% shiftq_yx = [-1.0 -0.5];
% padval = NaN;
% -----------------------------------------
% ____________________________________________________________
% 20.08.2018
% - floor+floor+round (with shift_yx=[0,0]: always ...,-1,0,1,... or ...,-2,-1,0,1,...)
% 03.07.2017
% - floor -> round
% 12.06.2017
% - shiftq_yx --> shift_yx (absolute shift in pixels)


[ny,nx] = size(mat);

if nargin<4; padval=0; end
if nargin<3 || isempty(shift_yx); shift_yx=[0 0]; end
if isempty(n_yx)
	noy = ny;
	nox = nx;
elseif length(n_yx)~=2
	error('Argument must be a list: n_yx=[noy,nox] ')
elseif isnan(n_yx(1)) && isnan(n_yx(2))
	error('Only one of n_yx=[noy, nox] can be NaN')
else
	if isnan(n_yx(1)); n_yx(1)=size(mat,1); end
	if isnan(n_yx(2)); n_yx(2)=size(mat,2); end
	noy = n_yx(1);
	nox = n_yx(2);
end

%y = round( ny/2 - noy/2 + (1:noy) + shift_yx(1) );
%x = round( nx/2 - nox/2 + (1:nox) + shift_yx(2) );
y = floor(ny/2) - floor(noy/2) + (1:noy) + round(shift_yx(1));
x = floor(nx/2) - floor(nox/2) + (1:nox) + round(shift_yx(2));

py1 = length(y(y<1));
py2 = length(y(y>ny));
px1 = length(x(x<1));
px2 = length(x(x>nx));

y=y(y>=1); y=y(y<=ny);
x=x(x>=1); x=x(x<=nx);

if length(y)==0; warning('No content left with shiftq_y=%d',round(shift_yx(1))); end
if length(x)==0; warning('No content left with shiftq_x=%d',round(shift_yx(2))); end

mat = mat(y,x);
mat = padarray(mat, [py1 0], padval, 'pre');
mat = padarray(mat, [py2 0], padval, 'post');
mat = padarray(mat, [0 px1], padval, 'pre');
mat = padarray(mat, [0 px2], padval, 'post');


