function [psnr] = psnr(sig,ref)

psnr = 10*log10(max(ref(:).^2)/sum(abs(sig(:)-ref(:)).^2));
