function [ind,BestModel] = GMM(data,max_components)

AIC = zeros(1,max_components);
GMModels = cell(1,max_components);

options = statset('MaxIter',1000);
rng(1); % For reproducibility

for k = 1:max_components
    GMModels{k} = fitgmdist(data,k,'Regularize',0.1,'Options',options);
    AIC(k)= GMModels{k}.AIC;
end

[minAIC,numComponents] = min(AIC);

BestModel = GMModels{numComponents};

P = posterior(BestModel,data);

[~,ind] = max(P');

end