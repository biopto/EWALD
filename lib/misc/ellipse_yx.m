function H = ellipse_yx(N1_N2, R1_R2, d1_d2)
% ELLIPSE_YX Create 2D ellipse of ones that can be used to limit
% information in the Fourier spectrum of a projection to the area of 
% the numerical aperture
% N1_N2   = [N1,N2] dimensions of the output array
% R1,R2   = [R1,R2] one half of the main axes of the ellipse
% d1,d2   = [d1,d2] shift vector (no wrapping)

N1 = N1_N2(1); N2 = N1_N2(2);
R1 = R1_R2(1); R2 = R1_R2(2);
if nargin<3
	d1=0; d2=0;
else
	d1 = d1_d2(1); d2 = d1_d2(2);
end

[NN1, NN2] = ndgrid(1:N1, 1:N2);
H = ((NN1-N1/2-1-d1)/R1).^2 + ((NN2-N2/2-1-d2)/R2).^2 <= 1.0;
