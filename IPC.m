function [fidxcluster,rmclusters,c,mean_clusters_mat,clusters]=IPC(data,alpha,plot_res)

tic;

if nargin < 3, plot_res=0; end

idxcluster{1} = clustering_corrH(data,alpha);

orig_idxcluster{1} = idxcluster{1};

clusters      = cell(length(idxcluster{1}),1);
rmclusters    = cell(length(idxcluster{1}),1);
for j = 1:length(idxcluster{1})
    clusters{j} = data(orig_idxcluster{1}{j},:);
end

mean_clusters_mat = cellfun(@mod_mean,clusters,'UniformOutput', false); 
mean_clusters_mat = cell2mat(mean_clusters_mat);

k=2;
con =1;
crit = 0;

while(con)
    
     idxcluster{k} = clustering_corrH(mean_clusters_mat,alpha);
    
     clusters      = cell(length(idxcluster{k}),1);
     for j = 1:length(idxcluster{k})
         orig_idxcluster{k}{j}  = horzcat(orig_idxcluster{k-1}{idxcluster{k}{j}});
         clusters{j} = data(orig_idxcluster{k}{j},:);
     end
    
    %% Prune
    for j = 1:length(clusters)
        PCorR = 1-pdist2(clusters{j},mod_mean(clusters{j}),'Spearman');
        if(any(PCorR<(alpha)))
            indrm = find(PCorR<alpha);
            rmclusters{k}{j}= clusters{j}(indrm,:);
            clusters{j}(indrm,:) = [];
            rmclustersidx{k}{j} = orig_idxcluster{k}{j}(indrm);
            orig_idxcluster{k}{j}(indrm) = [];
        end
    end
    
   mean_clusters_mat = cellfun(@mod_mean,clusters,'UniformOutput', false);
   mean_clusters_mat = cell2mat(mean_clusters_mat);
    
%    display(['prune',num2str(toc)]);
    
    %% Reassign
    
    if(~isempty(rmclusters{k}))
    for j = 1:length(rmclusters{k})
        if(~isempty(rmclusters{k}{j}))
        for p = 1:size(rmclusters{k}{j},1)
            PCAd = 1-pdist2(rmclusters{k}{j}(p,:),mean_clusters_mat,'Spearman');
            [max_PCorAd, ind_max] = max(PCAd);
            if(max_PCorAd<alpha)
                clusters{length(clusters)+1}  =  rmclusters{k}{j}(p,:);
                orig_idxcluster{k}{length(orig_idxcluster{k})+1} =  rmclustersidx{k}{j}(p);
            else
                clusters{ind_max}= [clusters{ind_max};rmclusters{k}{j}(p,:)];
                orig_idxcluster{k}{ind_max} = [orig_idxcluster{k}{ind_max},rmclustersidx{k}{j}(p)];
            end
        end
      end
    end
    end
   
    mean_clusters_mat = cellfun(@mod_mean,clusters,'UniformOutput', false);
    mean_clusters_mat = cell2mat(mean_clusters_mat);
    
%    display(['reassign',num2str(toc)]);
    
    %% Stopping Conditions
    
    %1. The indices of the clusters are identical
    %2. The number of iterations is >100
    %3. The number of jumping genes is less than 5%.
    
    s1 = isequal([orig_idxcluster{k}{:}], [orig_idxcluster{k-1}{:}]);
    
    for l = 1:length(orig_idxcluster{k})
    markov_ind(orig_idxcluster{k}{l}) = l;
    end

    for p = 1:length(orig_idxcluster{k-1})
    mcl_ind(orig_idxcluster{k-1}{p}) = p;
    end

    mem = zeros(size(data,1),size(data,1));
    mem_a = zeros(size(data,1),size(data,1));
    D = zeros(size(data,1),1);
    for i= 1:size(data,1)
    ind_g_i = find(markov_ind(i)==markov_ind);
    ind_g_i_2 = find(mcl_ind(i)==mcl_ind);
    mem(i,ind_g_i) = 1;
    mem_a(i,ind_g_i_2) = 1;
    D(i) = pdist2(mem(i,:),mem_a(i,:),'hamming');
    end
    
    c(k) = sum(1 - (D'*size(data,1))./(sum(mem)+sum(mem_a)))./size(data,1);
         
    if(k>3)
    
    m = mean(abs(diff(c((k-3):k))));
         
    if((k>100) || m<0.02)   
        con1 = 1;
         while(con1)
             
             k=k+1;
             
             idxcluster{k} = clustering_corrH(mean_clusters_mat,alpha);
             
             clusters      = cell(length(idxcluster{k}),1);
             for j = 1:length(idxcluster{k})
                 orig_idxcluster{k}{j}  = horzcat(orig_idxcluster{k-1}{idxcluster{k}{j}});
                 clusters{j} = data(orig_idxcluster{k}{j},:);
             end
             mean_clusters_mat = cellfun(@mod_mean,clusters,'UniformOutput', false);
             mean_clusters_mat = cell2mat(mean_clusters_mat);
             
             C = 1-pdist(mean_clusters_mat,'Spearman');
             
             con1 = max(max(C))>(alpha);
      
         end
         break;
    end
    
    end
     
     display(['No of cycles: ',num2str(k)]);

     k=k+1;
    
end

fidxcluster = orig_idxcluster{k};

% [s,ind] = sort(n_clusters,'descend');
%
% BC = tril(corrcoef(mean_clusters_mat'),-1);
%
% for b = 1:floor(length(fidxcluster)./30)
% h8=figure(b);
% for gen = 1:30
% subplot(5,6,gen)
% ind2 = ind(gen+(b-1).*30);
% plot(clusters{ind2}','-*b')
% hold on;
% plot(mean_clusters{ind2},'--r','LineWidth',1.5)
% xlim([0,size(clusters{ind2},2)])
% ylim([min(min(clusters{ind2}))-.05,max(max(clusters{ind2}))+.05])
% title([' Max Cor ',num2str(max((BC(:,ind2)))),'Cl ',num2str(ind2)])
% hold off;
% end
% print(h8, '-depsc', ['/Users/michellecarey/Dropbox/Cluster/Figures/H/',num2str(b),'ClusterH.eps']);
% end

end


