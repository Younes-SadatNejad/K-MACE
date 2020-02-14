function [cnc,IDX] = KMACE(data,mmin,mmax)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
% K-MACE is a method for determining the correct number of clusters (CNC)
% Reference Authors:
% Soosan Behshti (e-mail: soosan@ee.ryerson.ca), Edward Nidoy, Faizan Rahman, Cite:
% "K-MACE and Kernel K-MACE Clustering"
% IEEE ACCESS (2020)
% Website: https://www.ee.ryerson.ca/~soosan/
% Code Developement: Faizan Rahman , Younes Sadat-Nejad
% Contact Info: soosan@ee.ryerson.ca 
% Copy right February 2020
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
% ---- INPUT ---
%  data : Data input that you want to find the number of clusters for
%  mmin : Minimum number of clusters to experiment with (OPTIONAL : Default value is 1)
%  mmax : Maximum number of clusters to experiment with (OPTIONAL : Default value is 15)
% ---- OUTPUT ---
% cnc   : Numer of clusters

%   Example of Use:
%   [cnc]=KMACE(data,1,10) : Provides Number of clusters in data 
%   IDX  = 
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
if nargin < 2
    mmin=1;
    mmax=15;
end
if ~isempty(data)

%% data clustering
[N,dim] = size(data);

alphan = inline('nmj^(1/3)','nmj');
beta = N;

for m=mmin:mmax
    
    [IDX2(m,:)] = kmeans(data,m,'replicate',40,'EmptyAction','singleton');  % cluster the dataset using k-means
    
    for i=1:m
        aa2 = find(IDX2(m,:)==i);
        ni2(m,i) = length(aa2);
        est_center(i,:,m) = mean(data(aa2,:));
        cm2 = repmat(est_center(i,:,m),[ni2(m,i) 1]);
        
        temp_data = data(aa2,:);
        est_cov = cov(temp_data);
        ysmi2(m,i) = norm([cm2-temp_data],'fro').^2;    % eq 19        
        
        Ami = eye(ni2(m,i)) - (1/ni2(m,i))*ones(ni2(m,i));
        xxi(m,i) =  -(ni2(m,i))*(log(ni2(m,i)/N));
        
        aa2=[];
        cm2 =[];
        temp_data = [];
        
    end;
    ysm2(m,1) = (1/N)*sum(ysmi2(m,:));  % eq 20
    xx(m,1) = (1/(N))*sum(xxi(m,:),2);
end

%% calculating upperbound of ACE
for k=mmin:mmax
    cov_array = calc_cov(data,IDX2(k,:));   
    for  m=mmin:mmax
        for i=1:m
            aa2 = find(IDX2(m,:)==i);
            temp_var = cov_array(aa2,:);
            mw = ((ni2(m,i)-1)/ni2(m,i))*sum(sum(temp_var));
            alpha = alphan(ni2(m,i));
            varvar = 1*mltply(temp_var);
            
            bb= 2*sum(sum(varvar));
            aa = sum(sum(temp_var.^2,2));
            
            term2 = 2*(mw-ysmi2(m,i)) - (alpha^2)*(4/(dim*ni2(m,i))*sum(sum(temp_var,2)));  % term 2 of eq 61 
            term3 = (mw-ysmi2(m,i))^2 - ((alpha^2)/(ni2(m,i).^2))*((2*((ni2(m,i)-1).^2)*aa)+(bb));  % term 3 of eq 61
            
            if isnan(term2) || isnan(term3)
                upp_deltai2(m,i) = ysmi2(m,i);
                low_deltai2(m,i) = ysmi2(m,i);
                disp(['warning 1 problem at ' num2str(m) '- ' num2str(i) 'cluster ' num2str(ni2(m,i))]);
            else
                upp_deltai2(m,i) = (max(real(roots([1 term2 term3])))); % solving for eq 61
                low_deltai2(m,i) = (min(real(roots([1 term2 term3]))));
                if upp_deltai2(m,i)<0
                    upp_deltai2(m,i)=NaN;
                end
                if low_deltai2(m,i)<0
                    low_deltai2(m,i)=NaN;
                end
            end
            if ni2(m,i)==0
                der_zsmi2(m,i) = upp_deltai2(m,i) + sum(sum(temp_var));
                var_zsmi2(m,i) = 0;
                disp(['warning 2 problems at ' num2str(m) '- ' num2str(i) 'cluster ' num2str(ni2(m,i))]);
            else
                var_zsmi2(m,i) = ((2*((1/ni2(m,i))^2))* sum(sum(temp_var.^2,2))) + 2*((1/ni2(m,i))^2)*sum(sum(varvar)); % eq 14
                der_zsmi2(m,i) = upp_deltai2(m,i) + sum(sum(temp_var)).*(1/ni2(m,i));   % eq 13
            end
            
            varvar = []; % reseting this temporary array for the next calculation
            temp_var = [];
            
        end
        
        upp_delta(m,k) =(1/N)* sum(upp_deltai2(m,:),2);
        low_delta(m,k) =(1/N)* sum(low_deltai2(m,:),2);
        
        % now taking the total zsm
        der_zsm2(m) = (1/N)* sum(der_zsmi2(m,:));   % eq 17
        der_var_zsm2(m) = ((1/N)^2) * sum(var_zsmi2(m,:));  % eq 18
        
        der_upp_zsm2(m,k) = der_zsm2(m) + beta*sqrt(der_var_zsm2(m)) ;  % eq 34
        
        criterion(m,k) = der_upp_zsm2(m,k);
        
    end
end

%% calculating optimum cnc
zz=repmat([mmin:mmax],[mmax-mmin+1 1]);

criterion(1:mmin-1,:)=[];
criterion(:,1:mmin-1)=[];
[minval, minloc] = min(criterion,[],1); % find the minimum and its location for each zsm plot
diagonal = diag(criterion);


dd = [mmin:mmax];

holderr=[];
for i=1:length(minloc)-2
    if sum(isnan(criterion(minloc(i)+1:minloc(i)+2,i))~=0) %a minpoint is said to be invalid if there are NaN values beside it
       holderr = [holderr i];
    end
end

dd(holderr)=[];
minloc(holderr)=[];

lmq = find(abs(minloc-dd)==min(abs(minloc-dd)));

lmp = minloc(lmq);
ddd = dd(lmq);

for i=1:length(lmp)
    tt(i,:) = [criterion(lmp(i),ddd(i))' criterion(ddd(i),ddd(i))'];
    tt(i,:) = tt(i,:)/min(tt(i,:)) - 1;
    scc(i) = sum(tt(i,:));
end
[~,mmm] = min(scc);
cnc = lmp(mmm); % estimated cnc

IDX  = kmeans(data,cnc,'replicate',40,'EmptyAction','singleton');
center = est_center(1:cnc,:,cnc); %the cluster center
end