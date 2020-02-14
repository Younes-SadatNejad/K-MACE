function [kernel_matrix, sigma] = calc_sigma(data, sigma)

% Calculating the kernel distance matrix of the dataset data for a value of
% sigma

X = data;

kernel_matrix = zeros(size(X,1),size(X,1));

for r = 1:size(X,1)
    for s = 1:size(X,1)
        kernel_matrix(r,s) = ker_calc(X(r,:),X(s,:),sigma);
    end
end

end