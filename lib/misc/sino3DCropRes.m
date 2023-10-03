% Crop & resample a 3D sinogram (:,:,theta)
%
% Copyright: 2015, Piotr L. Makowski, p.makowski@mchtr.pw.edu.pl
%

function [SINO] = sino3DCropRes(SINO, qc, Qd)
% SINO      - projections stacked along 3rd dimension (:,:,proj)
% qc        - cropping factor (<=1)
% Qd        - downsampling factor (downsampling for >=1)

if qc>1; error('qc should be <=1'); end
if Qd<=0; error('Qd<=0 is forbidden'); end

if qc~=1; disp(['qc=' num2str(qc)]); end
if Qd~=1; disp(['Downsampling =' num2str(Qd)]); end

%%% crop
[n1, n2, n3] = size(SINO);
if ndims(SINO) == 3
    SINO = ndcrop(SINO, [n1*qc n2*qc n3]);
else
    % only 1 projection in the sinogram
    SINO = ndcrop(SINO, [n1*qc n2*qc]);
end

%%% resample
[n1, n2, n3] = size(SINO);
if Qd>1 % downsampling
    
    nn1 = round(n1/Qd);
    nn2 = round(n2/Qd);
    for th = 1:n3
        SINO(1:nn1,1:nn2,th) = imresize(squeeze(SINO(:,:,th)), [nn1 nn2], 'bicubic'); % bicubic: 4x4 averaging window 
    end
    SINO = SINO(1:nn1,1:nn2,:);
    
elseif Qd<1 % upsampling
    
%     SINO = ndresize(SINO, [n1/Qd n2/Qd n3], 'linear'); % cubic gives artifacts, imresize always better (bicubic best)
    
    nn1 = round(n1/Qd); p1 = ceil((nn1-n1)/2); % ceil: excessive space
    nn2 = round(n2/Qd); p2 = ceil((nn2-n2)/2);
    SINO = padarray(SINO, [p1, p2, 0]);
    for th = 1:n3
        SINO(1:nn1,1:nn2,th) = imresize(squeeze(SINO(p1+1:p1+n1, p2+1:p2+n2, th)), [nn1 nn2], 'bicubic'); % bilinear: dups higher freqs
    end    
    SINO = SINO(1:nn1,1:nn2,:); % cut excessive space
    
end


%%% even dimensions
[n1, n2, ~] = size(SINO);
if mod(n1,2); SINO(end,:,:)=[]; end
if mod(n2,2); SINO(:,end,:)=[]; end