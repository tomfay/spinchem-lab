function V = constructTwoSpinCoupling(A,n,m,g)
% constructs a two spin coupling operator
% A coupling tensor defining the coupling V = sum_ab S_na A_ab S_mb
% no assumptions about A are made
% n & m are indices of spins, assuming n<m
% g = (g_1, ... , g_N) is a vector of spin hilbert space dimensions for the
% N spin system

% total spin Hilbert space dimension
d = prod(g) ;
g_n = g(n) ;
g_m = g(m) ;
% dimension of subspaces excluding n & m
d_1 = prod(g(1:(n-1))) ;
d_2 = prod(g((n+1):(m-1))) ;
d_3 = prod(g((m+1):end)) ;

% empty V in lower dimensional subspace
V_trunc = sparse([],[],[],g(n)*d_2*g(m),g(n)*d_2*g(m));

% S_nx A_xb S_mb terms
sum_A_ab_S_mb = sparse([],[],[],g_m,g_m) ;
sum_A_ab_S_mb = sum_A_ab_S_mb + A(1,1) * spinOperatorX(g_m) ;
sum_A_ab_S_mb = sum_A_ab_S_mb + A(1,2) * spinOperatorY(g_m) ;
sum_A_ab_S_mb = sum_A_ab_S_mb + A(1,3) * spinOperatorZ(g_m) ;
V_trunc = kron(spinOperatorX(g_n),kron(speye(d_2),sum_A_ab_S_mb)) ;

% S_ny A_yb S_mb terms
sum_A_ab_S_mb = sparse([],[],[],g_m,g_m) ;
sum_A_ab_S_mb = sum_A_ab_S_mb + A(2,1) * spinOperatorX(g_m) ;
sum_A_ab_S_mb = sum_A_ab_S_mb + A(2,2) * spinOperatorY(g_m) ;
sum_A_ab_S_mb = sum_A_ab_S_mb + A(2,3) * spinOperatorZ(g_m) ;
V_trunc = V_trunc+ kron(spinOperatorY(g_n),kron(speye(d_2),sum_A_ab_S_mb)) ;

% S_nz A_zb S_mb terms
sum_A_ab_S_mb = sparse([],[],[],g_m,g_m) ;
sum_A_ab_S_mb = sum_A_ab_S_mb + A(3,1) * spinOperatorX(g_m) ;
sum_A_ab_S_mb = sum_A_ab_S_mb + A(3,2) * spinOperatorY(g_m) ;
sum_A_ab_S_mb = sum_A_ab_S_mb + A(3,3) * spinOperatorZ(g_m) ;
V_trunc = V_trunc+ kron(spinOperatorZ(g_n),kron(speye(d_2),sum_A_ab_S_mb)) ;

% construct the full V 
V = kron(speye(d_1),kron(V_trunc,speye(d_3))) ;

end