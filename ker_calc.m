function [result] = ker_calc(x1, x2, sigma)

% Calculating the distance between two points using the Gaussian kernel
% function

result = exp(-((norm(x1 - x2)).^2)/(2*sigma));
