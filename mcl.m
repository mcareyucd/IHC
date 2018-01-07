function[idxcluster] = mcl(data,delta)

%data = csvread('data1.csv',0,1); % get rid of the first column

D = 1-corr(data','type','Spearman');
N = size(data,1);
W = D < 1-delta;
W = W ./ repmat(sum(W),N,1); 
[xv xi] = max(W);
ms = xi; % ms the membership vector
tic
iter = 1;
while(true)
    disp(['iteration ' num2str(iter)])
    W = (W^2).^2;
    W = W ./ repmat(sum(W),N,1); 
    [xv xi] = max(W);
    if( sum(abs(ms-xi)) == 0) 
        break
    end
    ms = xi;
    iter = iter + 1;
end
toc

um = unique(ms);

for j = 1:length(um)
    ind = find(ms==um(j));
    ind2(ind)=j;
end

for j = 1:max(ind2)
    ind = find(ind2==j);
    idxcluster{j} = ind; 
end

end
