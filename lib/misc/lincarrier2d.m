function [S] = lincarrier2d(S, df_yx)
% [S] = lincarrier2d(S, df_yx)
%
%	Add linear phase to signal, so that fft2 is cirshifted by non-integer values
%
% %%%%% INPUT
% s      - complex 2d array
% df_yx  = [dfy, dfx] circshift values measured in frequency samples (fractional)
%
% -----------------------------------------------------------------------------
% %% EXAMPLE
% S = double(rgb2gray(imread('lena.png')));
% [Ny, Nx] = size(S);
% df_yx = [0 Nx/4];
% -----------------------------------------------------------------------------
%
% Copyright: 2017, Piotr L. Makowski, p.makowski@mchtr.pw.edu.pl
%


[Ny, Nx] = size(S);
y = ((1:Ny)-floor(Ny/2)-1)/(Ny/2);
x = ((1:Nx)-floor(Nx/2)-1)/(Nx/2);
[yy, xx] = ndgrid(y,x);

S = S .* exp(1i*pi*df_yx(1).*yy) .* exp(1i*pi*df_yx(2).*xx);
	