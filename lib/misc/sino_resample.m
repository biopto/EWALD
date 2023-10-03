% Copyright: 2019, Piotr L. Makowski, p.makowski@mchtr.pw.edu.pl

function sinostack3 = sino_resample(sinostack3, N, method)

if size(sinostack3,1)~=size(sinostack3,2)
	error('Sinogram must be rectangular')
end
N0 = size(sinostack3,1);
if N<N0
	for p=1:size(sinostack3,3)
		sinostack3(1:N,1:N,p) = imresize(sinostack3(:,:,p), [N N], method);
	end
	sinostack3 = sinostack3(1:N,1:N,:);
elseif N>N0
	sinostack3 = padarray(sinostack3(1:N0,1:N0,:), (N-N0)*[1 1 0], 0, 'post');
	for p=1:size(sinostack3,3)
		sinostack3(:,:,p) = imresize(sinostack3(1:N0,1:N0,p), [N N], method);
	end
end
