function H = electronHamiltonianQM(omega_1_ext , omega_2_ext, exchange, dipolar_tensor,  A_qm_1, A_qm_2, g_qm_1,g_qm_2) 

    % can be optimized massively by factoring out QM part, only needs to be
    % computed once.

    % generate the total fields experienced by both spins   
    omega_1 = omega_1_ext  ;
    omega_2 = omega_2_ext  ;
    
    % construct the hamiltonian
    H_SW = electronHamiltonian(omega_1 , omega_2, exchange, dipolar_tensor) ;
    
    % construct the S_1.A_1k.I_1k part of the hamiltonian
    s_x = spinOperatorX(2) ;
    s_y = spinOperatorY(2) ;
    s_z = spinOperatorZ(2) ;
    d_full = 4 * prod(g_qm_1) * prod(g_qm_2) ;
    d_1 = prod(g_qm_1) ;
    d_2 = prod(g_qm_2) ;
    N_qm_1 = length(g_qm_1) ;
    N_qm_2 = length(g_qm_2) ;
    
    H_1 = 1.0i*zeros(d_full) ;

    % This part needs to be optimised! Can factorise parts trivially to
    % improve.
       for k = 1:N_qm_1
          A = A_qm_1(1,1+3*(k-1)) * spinOperatorX(g_qm_1(k)) ...
              + A_qm_1(1,2+3*(k-1)) * spinOperatorY(g_qm_1(k)) ...
              + A_qm_1(1,3+3*(k-1)) * spinOperatorZ(g_qm_1(k)) ;
          A_full = kron(eye(prod(g_qm_1(1:(k-1)))),kron(A,eye(prod(g_qm_1((k+1):(N_qm_1))))) );
          H_1 = H_1 +kron(s_x,kron(eye(2), kron(A_full,eye(d_2))) );
          A = A_qm_1(2,1+3*(k-1)) * spinOperatorX(g_qm_1(k)) ...
              + A_qm_1(2,2+3*(k-1)) * spinOperatorY(g_qm_1(k)) ...
              + A_qm_1(2,3+3*(k-1)) * spinOperatorZ(g_qm_1(k)) ;
          A_full = kron(eye(prod(g_qm_1(1:(k-1)))),kron(A,eye(prod(g_qm_1((k+1):(N_qm_1))))) );
          H_1 = H_1 +kron(s_y,kron(eye(2), kron(A_full,eye(d_2))) );
          A = A_qm_1(3,1+3*(k-1)) * spinOperatorX(g_qm_1(k)) ...
              + A_qm_1(3,2+3*(k-1)) * spinOperatorY(g_qm_1(k)) ...
              + A_qm_1(3,3+3*(k-1)) * spinOperatorZ(g_qm_1(k)) ;
          A_full = kron(eye(prod(g_qm_1(1:(k-1)))),kron(A,eye(prod(g_qm_1((k+1):(N_qm_1))))) );
          H_1 = H_1 +kron(s_z,kron(eye(2), kron(A_full,eye(d_2))) );
       end

    
    % construct the S_2.A_2k.I_2k part of the hamiltonian
    H_2 = 1.0i*zeros(d_full) ;
    

       for k = 1:N_qm_2
          A = A_qm_2(1,1+3*(k-1)) * spinOperatorX(g_qm_2(k)) ...
              + A_qm_2(1,2+3*(k-1)) * spinOperatorY(g_qm_2(k)) ...
              + A_qm_2(1,3+3*(k-1)) * spinOperatorZ(g_qm_2(k)) ;
          A_full = kron(eye(prod(g_qm_2(1:(k-1)))),kron(A,eye(prod(g_qm_2((k+1):(N_qm_2))))) );
          H_2 = H_2 +kron(eye(2),kron(s_x, kron(eye(d_1),A_full)) );
          A = A_qm_2(2,1+3*(k-1)) * spinOperatorX(g_qm_2(k)) ...
              + A_qm_2(2,2+3*(k-1)) * spinOperatorY(g_qm_2(k)) ...
              + A_qm_2(2,3+3*(k-1)) * spinOperatorZ(g_qm_2(k)) ;
          A_full = kron(eye(prod(g_qm_2(1:(k-1)))),kron(A,eye(prod(g_qm_2((k+1):(N_qm_2))))) );
          H_2 = H_2 +kron(eye(2),kron(s_y, kron(eye(d_1),A_full)) );
          A = A_qm_2(3,1+3*(k-1)) * spinOperatorX(g_qm_2(k)) ...
              + A_qm_2(3,2+3*(k-1)) * spinOperatorY(g_qm_2(k)) ...
              + A_qm_2(3,3+3*(k-1)) * spinOperatorZ(g_qm_2(k)) ;
          A_full = kron(eye(prod(g_qm_2(1:(k-1)))),kron(A,eye(prod(g_qm_2((k+1):(N_qm_2))))) );
          H_2 = H_2 +kron(eye(2),kron(s_z, kron(eye(d_1),A_full)) );
       end

    
    % construct the full mixed QM-SW Hamiltonian
    H = kron(H_SW,eye(d_1*d_2)) + H_1 + H_2 ;
    
end