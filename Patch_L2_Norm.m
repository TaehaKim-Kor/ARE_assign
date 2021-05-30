%% Cost Function Set
%L2 Norm
function cost = Patch_L2_Norm(cmp,tgt)
cost = sum((cmp-tgt).^2,'all');
ele = size(tgt,1)*size(tgt,2);
cost = cost / ele;
end
