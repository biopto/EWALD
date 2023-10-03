% Match resolution of different data by resample & cropping/padding

function [MAT] = match_resolution(MAT, MAT_dx, new_dx, new_size)
% MAT       - 2D data
% MAT_dx    - sample size of MAT
% new_dx    - new sample size
% new_size  - new data matrix size

%%% resample
[n1, n2] = size(MAT);
Qd = MAT_dx/new_dx;

if Qd>1 % downsampling
    
    nn1 = round(n1/Qd);
    nn2 = round(n2/Qd);
    MAT(1:nn1,1:nn2) = imresize(MAT, [nn1 nn2], 'bicubic');
    MAT = MAT(1:nn1,1:nn2);
    
elseif Qd<1 % upsampling
    
    nn1 = round(n1/Qd); p1 = ceil((nn1-n1)/2); % ceil: excessive space
    nn2 = round(n2/Qd); p2 = ceil((nn2-n2)/2);
    MAT = padarray(MAT, [p1, p2]);
    MAT(1:nn1,1:nn2) = imresize(MAT(p1+1:p1+n1, p2+1:p2+n2), [nn1 nn2], 'bicubic');
    MAT = MAT(1:nn1,1:nn2); % cut excessive space

end

%%% crop
[n1, ~] = size(MAT);
if new_size<n1
    MAT = ndcrop(MAT, [new_size new_size]);
else
    MAT = padarray(MAT, [round((new_size-n1)/2) round((new_size-n1)/2)]);
end
%%% even dimensions
[n1, n2] = size(MAT);
if mod(n1,2); MAT(end,:)=[]; end
if mod(n2,2); MAT(:,end)=[]; end