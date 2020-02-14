function [output] = kkmeans_converge_mod5_1(data,k,K_mat)

% Kernel k-means algorithm

X = data;
rep = 1;
final_label = zeros(size(X,1),rep);
total_dist = zeros(rep,1);

for g = 1:rep
    
    dist = zeros(size(X,1),k);
    converged = 0;
    temp_K_mat = K_mat;
    temp_K_mat(temp_K_mat == 1) = 0;
    
    rows = size(X,1);
    columns = k;
    final_assignments = zeros(rows, columns);

    for u = 1:k
        for  v = 1:(size(X,1)/k)
            if v == 1
                [min_kdist, min_kdistind] = max(temp_K_mat(:));
                [kr,kc] = ind2sub(size(temp_K_mat),min_kdistind);
                temp_K_mat(kr,kc) = 0;
                temp_K_mat(kc,kr) = 0;
                if max(final_assignments(kr,:)) == 0
                    final_assignments(kr,u) = 1;
                end
                if max(final_assignments(kc,:)) == 0
                    final_assignments(kc,u) = 1;
                end
            else
                mindist_sum = zeros(size(X,1),1);
                for w = 1:size(X,1)
                    mindist_sum(w,1) = temp_K_mat(kr,w) + temp_K_mat(kc,w);
                end
                [temp_mindist_sum, temp_mindist_sumind] = max(mindist_sum);
                temp_K_mat(temp_mindist_sumind,:) = 0;
                temp_K_mat(:,temp_mindist_sumind) = 0;
                if max(final_assignments(temp_mindist_sumind,:)) == 0
                    final_assignments(temp_mindist_sumind,u) = 1;
                end    
            end
        end
        temp_K_mat(kr,:) = 0;
        temp_K_mat(:,kr) = 0;
        temp_K_mat(kc,:) = 0;
        temp_K_mat(:,kc) = 0;
    end
    
    
    Ck = sum(final_assignments,1);
    centroid = zeros(k,size(X,2));
    
    while ~converged
        
        
        for j = 1:k
            dist(:,j) = diag(K_mat) -(2/(Ck(j)))*sum(repmat(final_assignments(:,j)',size(X,1),1).*K_mat,2) + ...
                Ck(j)^(-2)*sum(sum((final_assignments(:,j)*final_assignments(:,j)').*K_mat));
        end
        oldfinal_assignments = final_assignments;
        final_assignments = (dist == repmat(min(dist,[],2),1,k));
        final_assignments = double(final_assignments);
        Ck = sum(final_assignments,1);
        
        for i = 1:k
            clust_idx = find(final_assignments(:,i) == 1);
            centroid(i,:) = mean(X(clust_idx,:));
        end
        
        if sum(sum(oldfinal_assignments~=final_assignments))==0
            converged = 1;
        end
        
    end
    
    f_label = zeros(size(X,1),1);
    
    for i = 1:size(final_assignments,1)
        f_label(i) = find(final_assignments(i,:) == 1);
    end
    
    emp_clust = zeros(k,1);
    dist_add = zeros(k,1);
    
    for i = 1:k
        emp_clust(i) = ismember(i,f_label);
        if emp_clust(i) == 0
            emp_idx = randperm(size(X,1));
            f_label(emp_idx(1)) = i;
        end
            temp1 = find(final_assignments(:,i) == 1);
            for u = 1:size(temp1,1)
                dist_add(i) = dist_add(i) + dist(temp1(u),i);
            end
        total_dist(g,1) = sum(dist_add);
    end
    
    final_label(:,g) = f_label;
   
end

[temp2,temp2ind] = min(total_dist);
output = final_label(:,temp2ind);

end

