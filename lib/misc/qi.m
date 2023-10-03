% Bovik, A.C., 2002. A universal image quality index. IEEE Signal Processing Letters, 9(3), pp.81�84.

function [qi] = qi(X1,X2) 

n = numel(X1);if numel(X2)~=n;error('Arguments must be of the same size.');end
m1 = mean(X1(:));
m2 = mean(X2(:));
s1 = sqrt( (1/(n-1)) * sum(sum(  (X1-m1).*(X1-m1) )) ) ;
s2 = sqrt( (1/(n-1)) * sum(sum(  (X2-m2).*(X2-m2) )) ) ;
s12 =    ( (1/(n-1)) * sum(sum(  (X1-m1).*(X2-m2) )) ) ;

qi = 4*s12*m1*m2/( (s1^2+s2^2)*(m1^2+m2^2) ) ;

end
