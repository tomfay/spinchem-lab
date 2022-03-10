function J = getIndexFromMultiIndex(multi_index,dims)
% transforms a total index for a multi-index for a tensor, where j_k =
% 1...d_k, into a single total index J = 1... prod_k d_k.
N = numel(dims) ;
J = multi_index(1) ;
for k = 2:N
    J = dims(k)*(J-1) + multi_index(k) ;
end
end