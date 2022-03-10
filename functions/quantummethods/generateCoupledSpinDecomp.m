function [g_blocks,block_weights] = generateCoupledSpinDecomp(g,N)

I = (g-1)/2 ;
K_max = N * I ;
N_blocks = floor(K_max)+1 ;
K_blocks = (1:N_blocks)+(K_max-N_blocks) ;
g_blocks = 2*K_blocks+1 ;
block_weights = zeros([1,N_blocks]) ;

if (mod(g,2)==1)
    n_I = ceil(g/2) ;
    block_weights(n_I) = 1 ;

    for n=2:N
        block_weights_new = block_weights ;
        for k = 1:((g-1)/2)
            block_weights_new = block_weights_new + [zeros([1,k]),block_weights(1:(end-k))] ;
            block_weights_new = block_weights_new + [block_weights((k+1):(end)),zeros([1,k])] ;
        end
        block_weights = block_weights_new ;
    end
else
    n_I = ceil(g/2) ;
    block_weights(n_I) = 1 ;
    for n = 2:N
        block_weights_new = zeros([1,N_blocks]) ;
        if mod(n,2)==0
            for k = 1:g 
                J = max([0,k-g/2]) ;
                K = max([0,g/2 - k]) ; 
                block_weights_new = block_weights_new + [zeros([1,J]),block_weights((1+K):(end-J)),zeros([1,K])] ;
            end
        else
            for k = 1:g 
                J = max([0,k-g/2-1]) ;
                K = max([0,g/2 - k+1]) ; 
                block_weights_new = block_weights_new + [zeros([1,J]),block_weights((1+K):(end-J)),zeros([1,K])] ;
            end
        end
        block_weights = block_weights_new ;
    end

end

end