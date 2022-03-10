function multi_index = getMultiIndexFromIndex(J,dims)
% transforms a total index for a multi-index tensor (indexed from 1 to
% d_tot = prod_k d_k) into a set of individual indices (indexed 1 to d_k)
    N = numel(dims) ;
    multi_index = zeros([1,N]) ;
    J_k = J ;
    for k = N:-1:2
        J_k_new = floor((J_k-1)/dims(k))+1 ;
        multi_index(k) = J_k - dims(k)*(J_k_new-1) ;
        J_k = J_k_new ;
    end
    multi_index(1) = J_k; 
end