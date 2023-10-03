% Copyright: 2018, Piotr L. Makowski, p.makowski@mchtr.pw.edu.pl

function [mat] = subarray_3d(mat, n_tab, shift_tab, padval)
% [mat] = subarray_3d(mat, n_tab, shift_tab, padval)
%
% Extract a shifted sub-array, fill unknown regions with padval
% Convention:
% 	even array/subarray: ...,-2,-1,0,1,...
% 	odd array/subarray:  ...,-1,0,1,...
%
% %%%%%%%%% INPUT
%	mat       - 3D array
%	n_tab     = [no1,no2,no3] output frame size;
%	            set NaN for one or two values to keep them;
%	            set n_tab=[] to keep all dimensions
% %%%%%%%%% <optional>
%	shift_tab = [D1,D2,D3] frame position shift with respect to
%	            input array center, in pixels
%	padval    - (default=0) padding value/method for padarray()
%
% Copyright: 2018, Piotr L. Makowski, p.makowski@mchtr.pw.edu.pl
%
% ----------------------------------------
% mat = cat(3, [1 2 3; 4 5 6], [1 2 3; 4 5 6]*10, [1 2 3; 4 5 6].*100)
% n_tab = [4 4 4];
% shift_tab = [0 0 1];
% padval = NaN;
% ----------------------------------------

[n1,n2,n3] = size(mat);

if nargin<4; padval=0; end
if nargin<3 || isempty(shift_tab); shift_tab=[0 0 0]; end

if isempty(n_tab)
	no1 = n1;
	no2 = n2;
elseif length(n_tab)~=3
	error('Argument must be a list: n_tab=[no1,no2,no3] ')
elseif sum(isnan(n_tab))==3
	error('Only one or two dimensions can be NaN')
else
	if isnan(n_tab(1)); n_tab(1)=size(mat,1); end
	if isnan(n_tab(2)); n_tab(2)=size(mat,2); end
	if isnan(n_tab(3)); n_tab(3)=size(mat,3); end
	no1 = n_tab(1);
	no2 = n_tab(2);
	no3 = n_tab(3);
end

% round!
i1 = floor(n1/2) - floor(no1/2) + (1:no1) + round(shift_tab(1));
i2 = floor(n2/2) - floor(no2/2) + (1:no2) + round(shift_tab(2));
i3 = floor(n3/2) - floor(no3/2) + (1:no3) + round(shift_tab(3));

p1a = length(i1(i1<1));
p1b = length(i1(i1>n1));
p2a = length(i2(i2<1));
p2b = length(i2(i2>n2));
p3a = length(i3(i3<1));
p3b = length(i3(i3>n3));

i1=i1(i1>=1); i1=i1(i1<=n1);
i2=i2(i2>=1); i2=i2(i2<=n2);
i3=i3(i3>=1); i3=i3(i3<=n3);

if isempty(i1); warning('No content left with shiftq_1=%d',round(shift_tab(1))); end
if isempty(i2); warning('No content left with shiftq_2=%d',round(shift_tab(2))); end
if isempty(i3); warning('No content left with shiftq_3=%d',round(shift_tab(3))); end

mat = mat(i1,i2,i3);
mat = padarray(mat, [p1a 0 0], padval, 'pre');
mat = padarray(mat, [p1b 0 0], padval, 'post');
mat = padarray(mat, [0 p2a 0], padval, 'pre');
mat = padarray(mat, [0 p2b 0], padval, 'post');
mat = padarray(mat, [0 0 p3a], padval, 'pre');
mat = padarray(mat, [0 0 p3b], padval, 'post');


