function [cnc,IDX] =Kernel_KMACE(data,mmax,min_sigma,max_sigma,incr)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
% Kernel_K-MACE is a method for determining the correct number of clusters (CNC)
% Reference Authors:
% Soosan Behshti, Edward Nidoy, Faizan Rahman,
% Cite:
% "K-MACE and Kernel K-MACE Clustering" IEEE ACCESS (2020)
% Website: https://www.ee.ryerson.ca/~soosan/
% Code Developement: Faizan Rahman Younes Sadat-Nejad
% Contact Info: soosan@ee.ryerson.ca 
% Copy right February 2020
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
% ---- INPUT ---
%  data : Data input that you want to find the number of clusters for
%  mmax : Maximum number of clusters to experiment with (Default value is 15)
%  min_sigma: Minimum range for sigma (Default value is -5)
%  max_sigma: Maximum range for sigma (Default value is 15)
%  incr: increment between sigma values (Default value is 1)
% ---- OUTPUT ---
% cnc   : Numer of clusters
% IDX   : Cluster labels for each sample in the dataset

%   Example of Use:
%   [cnc]=Kernel_KMACE(data,1,10,) : Provides Number of clusters in data 
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
if nargin < 2
    mmax=15;
    min_sigma=-5;
    max_sigma=15;
    incr=1;
end
if ~isempty(data)

%% k-MACE in feature space for all values of sigma
    
p = min_sigma;
t=(min_sigma+1):incr:(max_sigma-min_sigma);
mmin = 1;
for i = 1:length(t)-1
    sigma = 2^(p+i*incr);

    [C1,fin_idx(:,i),need] = kk_MACE(data(:,:,1),mmin,mmax,sigma);
    verge(:,i) = need;
end
    
%% Optimum Gaussian kernel parameter estimation
x = repmat([1:mmax]',[1 i]);
y = repmat([t],[mmax 1]);
[minval, minloc] = min(verge,[],1);

[wen, fen] = max(minval);
fang = minval(fen:end);
grad = diff(fang);
for i = 1:length(grad)-1
    if grad(i) > 0 && grad(i+1) > 0
        grad_sum(i) = 0;
    else
        grad_sum(i) = abs(grad(i)) + abs(grad(i+1));
    end
end
[~,maxir] = max(grad_sum);
maxir = maxir + fen;
cnc = minloc(maxir);
IDX = fin_idx(:,maxir);

end