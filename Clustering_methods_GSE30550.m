load('GSE30550.mat')

%Subject
sub = 1;

%Correlation Threshold
alpha = 0.10;

%IHC
[fidxcluster,rmclusters,c,mean_clusters_mat,clusters]=IHC(gexpst{sub},alpha);
% fidxcluster: index of the genes in each cluster
% rmclusters: clusters being removed at the prune step
% c: convergence %of genes that are consectiviely being assigned to the
% same cluster
% mean_clusters_mat: cluster centers
% clusters

%IPC
[fidxcluster,rmclusters,c,mean_clusters_mat,clusters]=IPC(gexpst{sub},alpha);
% fidxcluster: index of the genes in each cluster
% rmclusters: clusters being removed at the prune step
% c: convergence %of genes that are consectiviely being assigned to the
% same cluster
% mean_clusters_mat: cluster centers
% clusters

%MCL
[idxcluster] = mcl(gexpst{sub},alpha);
% idxcluster: index of the genes in each cluster

%Generalized Mixture Model Clustering
[BestModel,P] = GMM(gexpst{sub},20);

[m,ind] = max(P');

