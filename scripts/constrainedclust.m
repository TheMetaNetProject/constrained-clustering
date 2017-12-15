function i = GR_clusterer_BNCfeatures_GW(mode, betaScale)

basedir1 = '/u/metanet/clustering/constrained-clustering/data/'
addpath(basedir1)
suffix = 'BNC2000'
load([basedir1,'F-sims2000',suffix,'.mat'])
load([basedir1,'constraints',suffix,'.mat'])
constraints = constraints + 1; % because it's made in python, which has a zero index	
A = simmatrix;
%%%    A = exp(-dmat);
type = 'NJW';
nClusts=200;
iters = 50;
switch mode
case 1
    [clusts, distortion] = gcut2(A, nClusts, type,iters);
    filename = ['clusterslist_',type,suffix];
    save([basedir1,filename,'.mat'], 'clusts')
case 2
    tic()
%%%%    betaScale = 50;
    filename = ['SS_clusterslist_EN_BNCfeaturesGW_noun_WangDavidson',type,'_',int2str(betaScale),'e-2_',int2str(iters),'iterspos']
    Q = eye(size(A))*0.5;
%%%%%%    %Q(97,74)=1; Q(288,819)=1;Q(288,1397)=1; Q(51,847)=1; Q(68,161)=1; Q(68,1223)=1; Q(68,443)=1; Q(68,8)=1; Q(68,50)=1; Q(11,819)=1; Q(193,613)=1; Q(193,248)=1; Q(230,587)=1; Q(416,500)=1; Q(416,89)=1; Q(598,518)=1; Q(460,518)=1; Q(460,819)=1; Q(11,197)=1; Q(142,33)=1; Q(1124,1301)=1; Q(591,160)=1; Q(1510,866)=1; Q(1173,957)=1; Q(940,595)=1; Q(11,378)=1; Q(1124,245)=1; Q(308,1312)=1; Q(1280,662)=1; Q(288,1134)=1; Q(1099,791)=1; Q(30,791)=1;
%%%%%%    Q(97,139)=1; Q(11,288)=1; Q(11,460)=1; Q(11,51)=1; Q(288,460)=1; Q(288,51)=1; Q(460,51)=1; Q(416,591)=1; Q(598,11)=1; Q(1124,1127)=1; Q(1127,11)=1; Q(1510,1315)=1; Q(1514,1510)=1; Q(520,1000)=1; Q(280,11)=1; Q(1735,1280)=1; Q(845,507)=1; Q(967,507)=1; Q(97,36)=1; 
    for i = 1:size(constraints,1)
        Q(constraints(i,1),constraints(i,2)) = 0.5;
    end
    Q = Q + Q';
    [clusts,distortion] = wangdavidson(A,Q,betaScale/100, nClusts, 'shi', iters)
    toc()
    save([basedir1,filename,'.mat'], 'clusts')
case 3
    tic()
    beta = 1e1;
    filename = ['SS_clusterslist_EN_BNCfeaturesGW_noun_JiXuZhuNJW_1e1_',int2str(iters),'iters']
%%%%%%    constraints = [97 139; 11 288; 11 460; 11 51; 288 460; 288 51; 460 51; 460 51; 416 591; 598 11; 1124 1127; 1127 11; 1510 1315; 1514 1510; 520 1000; 280 11; 1735 1280; 845 507; 967 507; 97 36];
%%%%%%    constraints = [97 74; 288 819; 288 1397; 51 847; 68 161; 68 1223; 68 443; 68 8; 68 50; 11 819; 193 613; 193 248; 230 587; 416 500; 416 89; 598 518; 460 518; 460 819; 11 197; 142 33; 1124 1301; 591 160; 1510 866; 1173 957; 94 595; 11 378; 1124 245; 308 1312; 1280 662; 288 1134; 1099 791; 30 791];
    U = zeros(size(A,1), size(constraints,1));
    for i = 1:size(constraints,1)
        U(constraints(i,1),i) = 1;
        U(constraints(i,2),i) = -1;
    end
    [clusts,distortion] = jixuzhu(A,U,beta, nClusts, type, iters)
    toc()
    save([basedir1,filename,'.mat'], 'clusts')
end
noun = textscan(fopen([basedir1,'vocabBNC2000.txt']), '%s');
%load([basedir1,'normalized_GR_matrix_2000_BNCfeaturesGW.mat'])
noun = noun{1}(find(sum(A,2)>0),:);  % actually we need the original count matrix (not the similarity matrix) here ... (TODO)

clustlist = {};
for k = 1:length(clusts)
    vec = clusts{k};
    for j = 1:size(vec,2) 
        clustlist{k,j} = noun{vec(j)};
    end
end
idx = zeros(2000,1);
clustlist2 = {};
for k = 1:length(clusts)
    clustlist2{k,1}=['cluster',int2str(k),': '];
end
fid = fopen([basedir1,filename,'.txt'],'w');
for k = 1:length(clusts)
    vec = clusts{k};
    for j = 1:size(vec,2) 
        clustlist2{k,1} = [clustlist2{k,1},' ', noun{vec(j)}];
    end
    fwrite(fid,unicode2native([clustlist2{k,1}(1:end),10], 'UTF-8'),'uint8');
end
fclose(fid);

save([basedir1,filename,'.mat'], 'clustlist', 'distortion', 'clusts', 'clustlist2')