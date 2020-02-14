%% calculating the estimated ace based on a clustering result
function [est_cov] = calc_cov(data,idx)

[N dim] = size(data);
est_cov = zeros(N,dim);
for i=1:max(idx)
   aa = find(idx==i);
   ni = length(aa);
   
   if ni>1
       temp_data = data(aa,:);
       dd = sqrt(sum((temp_data - repmat(mean(temp_data),[ni,1])).^2,2));
       ddd = find(dd<3*median(dd));
       temp_data = temp_data(ddd,:);
       
       est_cov(aa,:) = repmat(eig(cov(temp_data))',[ni,1]);
   else
       est_cov(aa,:) = zeros(ni,dim);
   end
   temp_data=[];
   aa=[];
end
