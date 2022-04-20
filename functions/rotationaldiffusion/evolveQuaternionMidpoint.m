function q_dt = evolveQuaternionMidpoint(q_0,dt,N,k_B_T)

    w_1 = randn([3,1]) ;
    sqrt_N_1 = sqrtm(N(q_0)) ;
    omega = sqrt(4*k_B_T/dt) * sqrt_N_1 * w_1 ;
    q_dt = rotateQuaternionByOmega(q_0,omega,0.5*dt) ;
    w_2 = randn([3,1]) ;
    omega = sqrt(k_B_T/dt) * N(q_dt) * (sqrt_N_1 \ (w_1 + w_2)) ;
    q_dt = rotateQuaternionByOmega(q_0,omega,dt) ;

end