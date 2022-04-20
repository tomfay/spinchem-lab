function [q_ts,ts] = evolveQuaternionFreeTraj(q_0,n_steps,dt,N,k_B_T)
    ts = (0:n_steps)*dt ;
    q_ts = zeros([4,n_steps+1]) ;
    q_ts(:,1) = q_0 ;
    q_t = q_0 ;
    for k = 1:n_steps
        q_t = evolveQuaternionMidpoint(q_t,dt,N,k_B_T) ;
        q_ts(:,k+1) = q_t ;
    end
end