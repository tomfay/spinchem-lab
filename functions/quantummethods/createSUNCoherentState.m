function psi_coh = createSUNCoherentState(N)

    % creates an SU(N) coherent state
    % samples from delta(|c|^1 - 1) by sampling 2 N normally distributed
    % variables and projecting onto the surface with |c|^2 = 1 
    
    % sample the random variables
    psi_coh = randn([N,1]) + 1.0i * randn([N,1]) ;
    
    % normalise
    psi_coh = psi_coh * ( 1.0 /  norm(psi_coh) ) ;

end