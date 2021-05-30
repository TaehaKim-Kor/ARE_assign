%% Cost Function Set
%Pearson Correlation Coefficient
function cost = Patch_CC(cmp,tgt)
cost = [];
for i=1:size(cmp,3)
    restgt=reshape(tgt(:,:,i),[size(tgt(:,:,i),1)*size(tgt(:,:,i),2),1]);
    rescmp=reshape(cmp(:,:,i),[size(cmp(:,:,i),1)*size(cmp(:,:,i),2),1]);
    Cor=corrcoef(restgt,rescmp);
    ecost=Cor(2);
    cost = [cost,ecost];
end
cost=norm(cost);
end