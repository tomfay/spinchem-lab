function [psi_coh,Omegas] = createCoherentStateProduct(g)
    
    % number of spins to generate coherent state product for
    N = length(g) ;
    d = prod(g) ;
    % generate the angles, phis and thetas, sampled from surface of sphere
    [phis,thetas] = generateRandomAngles(N) ;
    Omegas = [ phis ; thetas ] ;
    
%     psi_coh = [1] ;
    psi_coh = (1.0+0.0i)*ones([d,1]) ;
    for  k = 1:N
        g_k = g(k) ;
        J_k = (g_k -1)/2 ;
%         j_ks = (1:g_k)' ;
        M_ks = (J_k:-1:(-J_k))' ;
        x = cos(thetas(k)*0.5) ;
        y = sin(thetas(k)*0.5)*exp(1.0i*phis(k)) ;
        psi_coh_k = 1.0i*zeros([g_k,1]) ;
        for j = 1:g_k 
            psi_coh_k(j) = sqrt(nchoosek(g_k-1,J_k+M_ks(j))) ...
                * (x^(J_k+M_ks(j))) * (y^(J_k-M_ks(j))) ;
        end
%         psi_coh = kron(psi_coh,psi_coh_k) ;
        d_l = prod(g(1:(k-1))) ;
        d_u = prod(g((k+1):end)) ;
        psi_coh_k_full = kron(kron(ones([d_l,1]),psi_coh_k),ones([d_u,1])) ;
        psi_coh = psi_coh.*psi_coh_k_full ;
    end

end