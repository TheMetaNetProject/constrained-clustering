basedir1 = '/u/metanet/clustering/en/constrainedclustering/'
addpath '/u/metanet/clustering/es/verbs/'
basedir2 =  '/u/metanet/clustering/en/GW/'
load([basedir2,'dmat_BNCfeaturesGW.mat'])

A = exp(-dmat);
[clusts, distortion] = gcut2(A, 200, 'shi');
clustlist = {};
save([basedir1,'clusterslist_EN_collapsed_noun_constrained.mat'])

noun = textscan(fopen('C:\Users\e4gutier\Dropbox\IARPA\Metaphor Extraction\nounlist.txt'), '%s');
load('normalized_GR_matrix_2000_BNCfeaturesGW.mat')
noun = noun{1}(find(sum(A,2)>0),:);

clustlist = {};
for k = 1:200 
    vec = clusts{k};
    for j = 1:size(vec,2) 
        clustlist{k,j} = verb(vec(j));
    end
end
for k = 1:200
    clustlist2{k,1}=['cluster',int2str(k),': '];
end
for k = 1:200
    vec = clusts{k};
    for j = 1:size(vec,2) 
        clustlist2{k,1} = [clustlist2{k,1},' ', verb{vec(j)}];
    end
end

save([basedir1,'clusterslist_ES_collapsed_verb.mat'], 'clustlist', 'distortion', 'clusts', 'clustlist2')