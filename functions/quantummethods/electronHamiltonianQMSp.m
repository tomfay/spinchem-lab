function H = electronHamiltonianQMSp(omega_1_ext , omega_2_ext, exchange, dipolar_tensor,  A_qm_1, A_qm_2, g_qm_1,g_qm_2) 

    % can be optimized massively by factoring out QM part, only needs to be
    % computed once.

    % generate the total fields experienced by both spins   
    omega_1 = omega_1_ext  ;
    omega_2 = omega_2_ext  ;
    
    % construct the hamiltonian Hamitlonian is + 2*exchange* S_1.S_2
    H_el = sparse(electronHamiltonian(omega_1 , omega_2, exchange, dipolar_tensor)) ;
    
    % construct the S_1.A_1k.I_1k part of the hamiltonian
    s_x = sparse(spinOperatorX(2)) ;
    s_y = sparse(spinOperatorY(2)) ;
    s_z = sparse(spinOperatorZ(2)) ;
    d_full = 4 * prod(g_qm_1) * prod(g_qm_2) ;
    d_1 = prod(g_qm_1) ;
    d_2 = prod(g_qm_2) ;
    N_qm_1 = length(g_qm_1) ;
    N_qm_2 = length(g_qm_2) ;
    
    H_1 = sparse([],[],[],d_full,d_full,0) ;

    % This part needs to be optimised! Can factorise parts trivially to
    % improve.
       for k = 1:N_qm_1
          A = sparse(A_qm_1(1,1+3*(k-1)) * spinOperatorX(g_qm_1(k)) ...
              + A_qm_1(1,2+3*(k-1)) * spinOperatorY(g_qm_1(k)) ...
              + A_qm_1(1,3+3*(k-1)) * spinOperatorZ(g_qm_1(k))) ;
          A_full = kron(speye(prod(g_qm_1(1:(k-1)))),kron(A,speye(prod(g_qm_1((k+1):(N_qm_1))))) );
          H_1 = H_1 +kron(s_x,kron(speye(2), kron(A_full,speye(d_2))) );
          A = sparse(A_qm_1(2,1+3*(k-1)) * spinOperatorX(g_qm_1(k)) ...
              + A_qm_1(2,2+3*(k-1)) * spinOperatorY(g_qm_1(k)) ...
              + A_qm_1(2,3+3*(k-1)) * spinOperatorZ(g_qm_1(k))) ;
          A_full = kron(speye(prod(g_qm_1(1:(k-1)))),kron(A,speye(prod(g_qm_1((k+1):(N_qm_1))))) );
          H_1 = H_1 +kron(s_y,kron(speye(2), kron(A_full,speye(d_2))) );
          A = sparse(A_qm_1(3,1+3*(k-1)) * spinOperatorX(g_qm_1(k)) ...
              + A_qm_1(3,2+3*(k-1)) * spinOperatorY(g_qm_1(k)) ...
              + A_qm_1(3,3+3*(k-1)) * spinOperatorZ(g_qm_1(k))) ;
          A_full = kron(speye(prod(g_qm_1(1:(k-1)))),kron(A,speye(prod(g_qm_1((k+1):(N_qm_1))))) );
          H_1 = H_1 +kron(s_z,kron(speye(2), kron(A_full,speye(d_2))) );
       end

    
    % construct the S_2.A_2k.I_2k part of the hamiltonian
    H_2 = sparse([],[],[],d_full,d_full,0) ;
    

       for k = 1:N_qm_2
%           A_qm_2(:,(1+3*(k-1)):(3+3*(k-1)))
          A = A_qm_2(1,1+3*(k-1)) * spinOperatorX(g_qm_2(k)) ...
              + A_qm_2(1,2+3*(k-1)) * spinOperatorY(g_qm_2(k)) ...
              + A_qm_2(1,3+3*(k-1)) * spinOperatorZ(g_qm_2(k)) ;
          A_full = kron(speye(prod(g_qm_2(1:(k-1)))),kron(A,speye(prod(g_qm_2((k+1):(N_qm_2))))) );
          H_2 = H_2 +kron(speye(2),kron(s_x, kron(speye(d_1),A_full)) );
          A = A_qm_2(2,1+3*(k-1)) * spinOperatorX(g_qm_2(k)) ...
              + A_qm_2(2,2+3*(k-1)) * spinOperatorY(g_qm_2(k)) ...
              + A_qm_2(2,3+3*(k-1)) * spinOperatorZ(g_qm_2(k)) ;
          A_full = kron(speye(prod(g_qm_2(1:(k-1)))),kron(A,speye(prod(g_qm_2((k+1):(N_qm_2))))) );
          H_2 = H_2 +kron(speye(2),kron(s_y, kron(speye(d_1),A_full)) );
          A = A_qm_2(3,1+3*(k-1)) * spinOperatorX(g_qm_2(k)) ...
              + A_qm_2(3,2+3*(k-1)) * spinOperatorY(g_qm_2(k)) ...
              + A_qm_2(3,3+3*(k-1)) * spinOperatorZ(g_qm_2(k)) ;
          A_full = kron(speye(prod(g_qm_2(1:(k-1)))),kron(A,speye(prod(g_qm_2((k+1):(N_qm_2))))) );
          H_2 = H_2 + kron(speye(2),kron(s_z, kron(speye(d_1),A_full)) );
       end

    
    % construct the full mixed QM-SW Hamiltonian
    H = kron(H_el,speye(d_1*d_2)) + H_1 + H_2 ;
    
end