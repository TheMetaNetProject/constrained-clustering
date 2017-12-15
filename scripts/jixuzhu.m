function [clusts, distortion] = jixuzhu(A,U,beta,nClusts, mode, Runs)
%
% Document Clustering with Prior Knowledge
%
% This function implements the constrainted spectral clustering algorithm
% of Ji, Xu, and Zhu (2006).
% Written by: E.D. Gutierrez (edg@icsi.berkeley.edu) using code by Lihi
% Zelnik-Manor and Yair Weiss
% Input:   
%       A = affinity matrix
%       U = constraint matrix where each column vector U encodes particular training pair: if V_i and V_j are in the same cluster, U(i,c)=1 and U(j,c)=-1 
%       beta = bias toward satisfying U - the greater beta is, the more the 
%              solution will be biased toward satisfying U.  
%       nClusts = k, the number of clusters
%       Runs = number of runs of k-means
%
%
%+%+%+%+%+
U = U';
N = size(A);
DD = 1./(sum(A)+eps);
if mode=='NJW' %For the Ng Jordan Weiss normalization
	DD = sqrt(DD);
	DD = diag(DD);
	L = eye(N) - DD*(A - beta*U'*U)*DD;
else % For the Shi et al. normalization
	DD = diag(DD)
	L = eye(N) - DD*(A - beta*U'*U);
end
%%%%%%% Compute eigenvectors

    opts.issym = 1;
    opts.isreal = 1;
    opts.disp = 0;
    [V,~] = eigs(L,nClusts,0,opts); %COMPUTING THE SMALLEST EIGENVALUES!!!
    V = V(:,1:nClusts);    

%%%%%%% Normalize rows of V - NOT SURE IF THIS SHOULD BE DONE FOR JI ET AL. (PROBABLY)
for i=1:size(V,1);
    V(i,:)=V(i,:)/(norm(V(i,:))+1e-10);
end

%%%%%%%%%%%%%%%%%% Kmeans
%%%%%% Try Runs runs of k-means and save the one with minimal distortion
bestC = {};           %% a variable to keep the best clustering so far
bestD = 1/eps;        %% a variable to remember the best distortion so far
for nRuns=1:Runs
    %%%%%% Initialize centers 
    %% First center is set to one entry picked randomly 
    [~,pp] = max(rand(size(V,1),1)); 
    mu(1,:) = V(pp,:);
    %% The other centers are selected to be farthest from previous centers
    for i=2:nClusts
        ip = V*mu';
        minip = max(abs(ip'));
        [~,ii] = min(minip);
        mu(i,:) = V(ii,:);
    end
    %%%%%%%%%%% and now run K means
    for tt=1:15   %% 10 iterations for kmeans
        distM = dist2(V,mu);        %% initialize distance between points and centers
        [~,ii] = min(distM');      %% assign points to nearest center

        distort = 0;
        distort_across = 0;
        clear clusts;
        for nn=1:nClusts
            I = find(ii==nn);       %% indices of points in cluster nn
            J = find(ii~=nn);       %% indices of points not in cluster nn
            clusts{nn} = I;         %% save into clusts cell array
            if (length(I)>0)
                mu(nn,:) = mean(V(I,:));               %% update mean
                %% Compute within class distortion
                muB = repmat(mu(nn,:),length(I),1);
                distort = distort+sum(sum((V(I,:)-muB).^2));
                %% Compute across class distortion
                muB = repmat(mu(nn,:),length(J),1);
                distort_across = distort_across + sum(sum((V(J,:)-muB).^2));
            end
        end
        %% Set distortion as the ratio between the within
        %% class scatter and the across class scatter
        distort = distort/(distort_across+eps);
        if (distort<bestD)   %% save result if better than the best so far
            bestD=distort;
            bestC=clusts;
        end
    end
end

%% Finally, delete empty clusters
pp=1;
for nn=1:nClusts
    if (length(bestC{nn})>0)
        clusts{pp} = bestC{nn}; %#ok<AGROW>
        pp = pp+1;
    end
end
distortion  = bestD;