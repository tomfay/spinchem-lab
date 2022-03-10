function A = sumSpinOpsX(g ,a )
    % number of spins
    N = length(g) ;
    % total dimension of spin hilbert space
    d = prod(g) ;
    % empty array for A = sum_k a_k I_kx
    A = sparse([],[],[],d,d,0) ;
    for k = 1:N
        d_1 = prod(g(1:(k-1))) ;
        d_2 = prod(g((k+1):(N))) ;
        A = A + kron(speye(d_1),kron(a(k)*spinOperatorX(g(k)),speye(d_2))) ; 
    end

end