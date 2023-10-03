% Plots: extract smaller cube from REC
function [M, c1,c2,c3] = center_cube_3d(M, c1,c2,c3, Nc, padv)
% M			- cube array
% c1,c2,c3  - index coords of the ceter point
% Nc		- size of the output cube centered around (c1,c2,c3)
% padv		- padding value for parts outside M

if nargin<6; padv=0; end
if ~isempty(M) && ~isnan(Nc)
	[n1, n2, n3] = size(M);
	M = subarray_3d(M, [Nc Nc Nc], ([c1 c2 c3]-[n1 n2 n3]/2), padv);
	c1 = Nc/2+1;
	c2 = Nc/2+1;
	c3 = Nc/2+1;
end