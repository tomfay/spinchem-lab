function psi_new = propagateSIA(H,psi_0,dt,p) 
    % approximates exp(-i H dt) psi_0 using the short iterative arnoldi
    % method.
    % p is the krylov subspace dimension
    
    % generate the kyrlov subspace representation of H
    [H_krylov, krylov_basis] = generateKrylovSubspace(H, psi_0, p) ;
    
    % construct the new state
    psi_0_krylov = zeros(p,1) ;
    psi_0_krylov(1) = norm(psi_0) ;
    psi_new = krylov_basis * ( expm(-1.0i * dt * H_krylov) *  psi_0_krylov) ;
    
    
end