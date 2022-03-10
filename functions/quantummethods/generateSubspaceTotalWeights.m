function [total_weights,total_weight_Zs, Z_subspaces] = generateSubspaceTotalWeights(block_weights,n_subspaces,g_blocks)
    N_subspaces_tot = prod(n_subspaces) ;
    total_weights = ones([1,N_subspaces_tot]) ;
    total_weight_Zs = ones([1,N_subspaces_tot]) ;
    Z_subspaces = ones([1,N_subspaces_tot]) ; 
    n_blocks = numel(n_subspaces) ;

    for J = 1:N_subspaces_tot
        j = getMultiIndexFromIndex(J,n_subspaces) ;
        Z_subspace = 1 ;
        for k = 1:n_blocks
            total_weights(J) = total_weights(J) * block_weights{k}(j(k)) ;
            Z_subspace = Z_subspace *  g_blocks{k}(j(k)) ;
        end
        total_weight_Zs(J) = total_weights(J) * Z_subspace ;
        Z_subspaces(J) = Z_subspace ;
    end

end