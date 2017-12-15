function [clusts, distortion] = wangdavidson(A,Q,betaScale,nClusts, mode, Runs)
%
% Flexible Constrained Spectral Clustering
%
% This function implements the constrainted spectral clustering algorithm
% of Wang and Davidson (2010).
% Written by: E.D. Gutierrez (edg@icsi.berkeley.edu) using code by Lihi
% Zelnik-Manor and Yair Weiss
% Input:   
%       A = affinity matrix
%       Q = constraint matrix
%       beta = bias toward satisfying Q - the greater beta is, the more the 
%              solution will be biased toward satisfying Q.  Set beta 
%              between lambda_min*vol_G and lambda_max*vol_G
%       nClusts = k, the number of clusters
%
%
%
%+%+%+%+%+


c = 10; %rule of thumb of how many times more eigenvectors we need to ...
       %compute to find the k largest ones with positive eigenvalues
vol_G = sum(sum(A));
N = size(A,1);
DD = 1./(sum(A)+eps);
if strcmp(mode,'NJW')
      'NJW'
	DD = sqrt(DD);
	DD = diag(DD);
	Lbar = eye(N) - DD*A*DD;
	Qbar = DD*Q*DD;
%	Qbar = Q;

else % For the Shi et al. algorithm
	DD = diag(DD);
	Lbar = eye(N) - DD*A;
	Qbar = DD*Q;
end
lambda_max = eigs(Qbar, 1); %MAKE SURE THIS IS COMPUTING THE LARGEST EIGVAL OF QBAR
try 
    lambda_min = eigs(Qbar, 1, 'sm');
catch
    lambda_min = 0;
end

beta = (betaScale*lambda_min+(1-betaScale)*lambda_max)*vol_G  %Made-up rule-of-thumb
%beta = betaScale*lambda_max *vol_G *(0.5+0.4*sum(sum(Q~=0))/(size(Q,1)*size(Q,2)))

if beta >= lambda_max*vol_G
    u_star = [];
else
    k = min(size(Lbar,2)-2, c*nClusts);
%    Qbar = eye(size(A));    'Warning this is code with manipulated qbar'
    [vecs, vals] = eigs(Lbar, (Qbar - beta/vol_G*eye(N)), k, 'lm');  % MAKE SURE THE SIGN OF Lbar IS RIGHT HERE
    vecs;
    nClusts
    sizevecs = size(vecs)
    V = zeros(N, c*nClusts);
    i = 1;
    j = 1;
    while i<nClusts
        if j>=size(Lbar,2)-3
                 i
  		 i=1e5;
              
        end
        j = j + 1;

        if vals(j,j)>0
            i = i+1;
            V(:,i) = vecs(:, j)*vol_G/norm(vecs(:,j)+1e-10);
        end
    end
    
end
for i=1:size(V,2)
    V(:,i) = DD*V(:,i);
end
vals(j,j)
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
clear clusts
for nn=1:nClusts
    if (length(bestC{nn})>0)
        clusts{pp} = bestC{nn}; %#ok<AGROW>
        pp = pp+1;
    end
end
distortion  = bestD;