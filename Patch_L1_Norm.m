%% Cost Function Set
%L1 Norm
function cost = Patch_L1_Norm(cmp,tgt)
cost = sum(abs(cmp-tgt),'all');
ele = size(tgt,1)*size(tgt,2);
cost = cost / ele;
end