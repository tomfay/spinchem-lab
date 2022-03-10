function [block_weights,g_blocks,n_subspaces] = createSymmetryBlocks(g,N)
    n_blocks = numel(N) ;
    block_weights = cell([1,n_blocks]) ;
    g_blocks = cell([1,n_blocks]) ;
    n_subspaces = zeros([1,n_blocks]) ;
    % generate the set of symmetry block weights and subspace dimensions
    % for each symmetry block
    for k = 1:n_blocks
        [g_blocks_k,block_weights_k] = generateCoupledSpinDecomp(g(k),N(k)) ;
        g_blocks{k} = g_blocks_k ;
        block_weights{k} = block_weights_k ;
        n_subspaces(k) = numel(g_blocks_k) ;
    end


end