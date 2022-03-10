function H = electronHamiltonianQMGen(omega_1_ext , omega_2_ext, exchange, dipolar_tensor,  A_qm_1, A_qm_2, g_qm_1,g_qm_2) 

    % can be optimized massively by factoring out QM part, only needs to be
    % computed once.

    % generate the total fields experienced by both spins   
    omega_1 = omega_1_ext  ;
    omega_2 = omega_2_ext  ;
    
    % construct the hamiltonian Hamitlonian is + 2*exchange* S_1.S_2
    H_el = sparse(electronHamiltonian(omega_1 , omega_2, exchange, dipolar_tensor)) ;
    
    % some useful quantities
    g = [2,2,g_qm_1,g_qm_2] ;
    g_nuc = prod([g_qm_1,g_qm_2]) ;
    N_1 = length(g_qm_1) ;
    N_2 = length(g_qm_2) ;
    
    % full Hamiltonian with just the electron spin term
    H = kron(H_el,speye(g_nuc)) ;

    % spin indices for the radical 1 spins
    m_list = (2+1):(2+N_1) ;
    for k = 1:N_1
       A = A_qm_1(:,(3*(k-1)+1):(3*k))  ;
       n = 1 ;
       m = m_list(k) ;
       H = H + constructTwoSpinCoupling(A,n,m,g) ;
    end
    
     % spin indices for the radical 1 spins
    m_list = (2+N_1+1):(2+N_1+N_2) ;
    for k = 1:N_2
       A = A_qm_2(:,(3*(k-1)+1):(3*k))  ;
       n = 2 ;
       m = m_list(k) ;
       H = H + constructTwoSpinCoupling(A,n,m,g) ;
    end
    
    
end