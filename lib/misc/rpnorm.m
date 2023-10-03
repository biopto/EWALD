% relative p-norm
function [e] = rpnorm(s,s0,p)
% s  - signal
% s0 - reference
% p  - norm number (1=MAE, 2=RMSE)

e = norm(s(:)-s0(:),p)/norm(s0(:),p);
