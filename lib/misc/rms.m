function [r]=rms(X)

r = sqrt(sum(abs(X).^2)/(numel(X)-0));
