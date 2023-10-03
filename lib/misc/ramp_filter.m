% Strict ramp filter implementation accounting for finite pixel size
% (gridrec paper, Marone & Stampanoni, 2012)
%
% Copyright: 2018, Piotr L. Makowski, p.makowski@mchtr.pw.edu.pl
%

function [R] = ramp_filter(fs, mode, shift)
% fs	- vector of frequencies to filter (even, left(-) arm longer)
% mode	- (=1, default): 1D ramp filter, (=2): 2D conical filter
% shift = [d1,d2] (default=[0 0]) offset in units of ramp pixel
% R		- 1D or 2D array, depending on mode
if nargin<3; shift = [0 0]; end
if nargin<2; mode=1; end

% N = 100;
N = length(fs);
sigma = zeros(1,N);
T = (1:N)-floor(N/2)-1;
ids = 2-mod(floor(N/2),2) : 2 : N; % odd pixels counting from "0"=floor(N/2)+1
sigma(ids) = -1./(pi*T(ids)).^2;
sigma(floor(N/2)+1) = 1/4;


% 1d
R = fftshift(fft(ifftshift(sigma)));
R = R/(abs(R(end))); % correct with even N
R = R(:);

% % ----------------------
% ids = 1:length(f);
% % ids = find(f>=-5);
% % ids = find(f>=0 & f<1e-3);
% fi = f(ids);
% Ri = R(ids);
% 
% figure(1)
% plot(fi, abs(Ri),'o-b'); hold on
% plot(fi, imag(Ri),'xr'); hold off
% legend({'abs','imag'})
% % xlim([0 1
% % ylim([0 1e-3])
% % ----------------------


% 2d
if mode==2
    f = ((1:N)-floor(N/2)-1)*1;
    Rfun = @(ro) interp1(f,R,ro,'linear',1);
%     % ----------------------
%     figure(11)
%     plot(f, abs(R),'o-b'); hold on
%     plot(interp(f,10),Rfun(interp(f,10)),'-g'); hold off
%     legend({'points','interpolation'})
%     % ----------------------
    
    [fx,fy] = meshgrid(f,f);
    R = Rfun(sqrt(fx.^2+fy.^2));
	
	if(sum(abs(shift))>0)
		fR = fftshift(ifft2(ifftshift(R)));
% 			figure;imagesc(log10(1e-2+abs(fR))); axis image; colormap jet;
		fR = lincarrier2d(fR, shift);
		R = fftshift(fft2(ifftshift(fR)));
% 			figure; imagesc(real(R)); axis image; colormap jet; colorbar; title('real(R)')
% 			figure; plot(real(R(N/2+round(-shift(1)),:)))
	end
% 		figure; mesh(fx,fy,real(R));
end

R = real(R);
