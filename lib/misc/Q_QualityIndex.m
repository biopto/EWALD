% Algorithm implemented from:
% Bovik, A.C., 2002. A universal image quality index. IEEE Signal Processing Letters, 9(3), pp.81�84.
% Author: Wojciech Krauze

function [Q] = Q_QualityIndex(X,Y) 
% mean_X = mean(mean(X)) ;
% mean_Y = mean(mean(Y)) ;

N = size(X,1)*size(X,2) ;

mean_X = sum(X(:))/N ;
mean_Y = sum(Y(:))/N ;

sigmax = sqrt( (1/(N-1)) * sum(sum(  (X-mean_X).*(X-mean_X) )) ) ;
sigmay = sqrt( (1/(N-1)) * sum(sum(  (Y-mean_Y).*(Y-mean_Y) )) ) ;
sigmaxy =    ( (1/(N-1)) * sum(sum(  (X-mean_X).*(Y-mean_Y) )) ) ;

Q = 4*sigmaxy*mean_X*mean_Y/( (sigmax^2+sigmay^2)*(mean_X^2+mean_Y^2) ) ;

end