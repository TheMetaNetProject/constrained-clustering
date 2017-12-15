function [clusts, distortion, k] = xianggong(A,kMax, mode1, EMIters)
%%+%+%+%+%+%+%+%+%+%+%+%+%+%+%+%+%+%+%+%+%+%+%+%+%+%+%+%+%+%+%+%+%+%+%+%+%+
% Spectral Clustering with Eigenvector Selection
%
% This function implements the constrainted spectral clustering algorithm
% of Xiang and Gong (2007).
% Written by: E.D. Gutierrez (edg@icsi.berkeley.edu) using code by Lihi
% Zelnik-Manor and Yair Weiss
% Input:   
%       A = affinity matrix
%       kMax = maximum size of k, the number of clusters
%
%
%
%+%+%+%+%+%+%+%+%+%+%+%+%+%+%+%+%+%+%+%+%+%+%+%+%+%+%+%+%+%+%+%+%+%+%+%+%+%
%+%+%+%+%+%+%+%+%+%+%+%+%+%+%+%+%+%+%+%+%+%+%+%+%+%+%+%+%+%+%+%+%+%+%+%+%+%

N = size(A,1);
DD = 1./(sum(A)+eps);
if mode1=='NJW' % For Ng Jordan Weiss normalization; (used in orig. paper)
	DD = sqrt(DD);
	DD = diag(DD);
	L = eye(N) - DD*A*DD;
else % For the Shi et al. normalization
	DD = diag(DD);
	L = eye(N) - DD*A;
end%%%%%%% Compute eigenvectors
if (useSparse)
    opts.issym = 1;
    opts.isreal = 1;
    opts.disp = 0;
    [V,~] = eigs(L,kMax,'LM',opts); %COMPUTING THE LARGEST EIGENVALUES!!!
%     [VV,ss]=svds(L,kMax,1,opts);
else
    opts.issym = 1;
    opts.isreal = 1;
    opts.disp = 0;
    [V,~] = eigs(L,kMax,'LM',opts);
    %[V,~] = svd(L); %MAKe sure THIS IS COMPUTING THE LARGEST EIGENVALUES!! 
end

%%%%%%%%%%%%%%%%%%%Eigenvector Relevance Computation%%%%%%%%%%%%%%%%%%%%%%%
VV = [];
for k = 1:kMax
    e_k = V(:,k);
    mu_k1 = (1/N)*sum(e_k);
    sigma_k1 = (1/N)*sum((e_k-mu_k1).^2);
    %%%%INITIALIZATION%%%%
    R_e_k = 0.5;
    w_k = rand(1);
    mu_k2 = rand(1);    mu_k3 = rand(1);
    sigma_k2 = rand(1);    sigma_k3 = rand(1);
    for j = 1:EMIters
    %%%%%%%%E-STEP%%%%%%%%
        for n = 1:N
            u = (1-R_e_k)*normpdf(e_k(n), mu_k1, sigma_k1);
            v = w_k*R_e_k*normpdf(e_k(n), mu_k2, sigma_k2);
            w = (1-w_k)*R_e_k*normpdf(e_k(n), mu_k3, sigma_k3);
		    h1_k(n) = u/(u+v+w+eps);
		    h2_k(n) = v/(u+v+w+eps);
		    h3_k(n) = w/(u+v+w+eps);
        end
   %%%%%%%%M-STEP%%%%%%%%
        R_e_k = 1 - (1/N)*sum(h1_k);
        w_k = 1/((R_e_k*N)*sum(h2_k)+eps);
        mu_k2 =sum(h2_k.*e_k')/sum(h2_k+eps);mu_k3 = sum(h3_k.*e_k')/sum(h3_k+eps);
        sigma_k2 = sum(h2_k.*(e_k'-mu_k2).^2)/sum(h2_k+eps);  sigma_k3 = sum(h3_k.*(e_k'-mu_k3).^2)/sum(h3_k+eps);
    end
    if R_e_k > 0.5
        VV = [VV, R_e_k*e_k];
    end
end


for k = 2:kMax
    switcher = 0;
    while switcher<200
        [k, switcher]
        clear obj
        try
            obj = gmdistribution.fit(VV,k, 'Options', statset('MaxIter',400));
            switcher = switcher + 1;
            if switcher ==1
                AIC(k) = obj.AIC;
                BIC(k) = obj.BIC;
            end
            if AIC(k) > obj.AIC
                AIC(k) = obj.AIC;
            end
            if BIC(k) > obj.BIC
                BIC(k) = obj.BIC;
            end
        catch
        end
    end
end

k = find(BIC==min(BIC))
clusts = 0;
distortion = 0;
%%%%%%% Normalize rows of V - NOT SURE IF THIS SHOULD BE DONE FOR JI ET AL. (PROBABLY)
% for i=1:size(V,1);
%     V(i,:)=V(i,:)/(norm(V(i,:))+1e-10);
% end
% 
% %%%%%%%%%%%%%%%%%% Kmeans
% %%%%%% Try 500 runs of k-means and save the one with minimal distortion
% bestC = {};           %% a variable to keep the best clustering so far
% bestD = 1/eps;        %% a variable to remember the best distortion so far
% for nRuns=1:500
%     %%%%%% Initialize centers 
%     %% First center is set to one entry picked randomly 
%     [~,pp] = max(rand(size(V,1),1)); 
%     mu(1,:) = V(pp,:);
%     %% The other centers are selected to be farthest from previous centers
%     for i=2:kMax
%         ip = V*mu';
%         minip = max(abs(ip'));
%         [~,ii] = min(minip);
%         mu(i,:) = V(ii,:);
%     end
%     %%%%%%%%%%% and now run K means
%     for tt=1:15   %% 10 iterations for kmeans
%         distM = dist2(V,mu);        %% initialize distance between points and centers
%         [~,ii] = min(distM');      %% assign points to nearest center
% 
%         distort = 0;
%         distort_across = 0;
%         clear clusts;
%         for nn=1:kMax
%             I = find(ii==nn);       %% indices of points in cluster nn
%             J = find(ii~=nn);       %% indices of points not in cluster nn
%             clusts{nn} = I;         %% save into clusts cell array
%             if (length(I)>0)
%                 mu(nn,:) = mean(V(I,:));               %% update mean
%                 %% Compute within class distortion
%                 muB = repmat(mu(nn,:),length(I),1);
%                 distort = distort+sum(sum((V(I,:)-muB).^2));
%                 %% Compute across class distortion
%                 muB = repmat(mu(nn,:),length(J),1);
%                 distort_across = distort_across + sum(sum((V(J,:)-muB).^2));
%             end
%         end
%         %% Set distortion as the ratio between the within
%         %% class scatter and the across class scatter
%         distort = distort/(distort_across+eps);
%         if (distort<bestD)   %% save result if better than the best so far
%             bestD=distort;
%             bestC=clusts;
%         end
%     end
% end

%% Finally, delete empty clusters
% pp=1;
% for nn=1:kMax
%     if (length(bestC{nn})>0)
%         clusts{pp} = bestC{nn}; %#ok<AGROW>
%         pp = pp+1;
%     end
% end
% distortion  = bestD;