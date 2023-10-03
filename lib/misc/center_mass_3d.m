% Plots: find center of mass in REC
function [c1, c2, c3] = center_mass_3d(M)
% M	        - 3d image with zero background
% c1,c2,c3  - index coords of center of mass of abs(M)

med=median(M(:));
M=M-med;
M(M<0)=0;

S1 = sum(sum(abs(M),3),2);
S2 = sum(sum(abs(M),3),1);
S3 = squeeze(sum(sum(abs(M),2),1));
	S3=S3.*tukeywin(length(S3),0.25);% supress edges
c1 = find(S1==max(S1));
c2 = find(S2==max(S2));
c3 = mean(find(S3==max(S3)));% in case of very flat profile