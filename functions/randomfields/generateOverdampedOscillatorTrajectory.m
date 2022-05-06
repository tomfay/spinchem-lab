function q = generateOverdampedOscillatorTrajectory(q_0, alpha, k, k_B_T, dt, n_time)

    % equation of motion is
    % dq(t)/dt = -(k/alpha) q(t) + (1/alpha) beta(t)
    % beta is a random force, alpha is the friction coefficent, k is the force constant.

    % calculate k_n * dt / 2 alpha_n 
    k_dt_over_2_alpha = (k ./ alpha) * dt*0.5 ;
    sigma_beta = sqrt(2.0*k_B_T * dt./alpha) ;

    q = zeros([size(q_0,1),n_time+1]) ;
    q(:,1) = q_0 ;
    q_t = q_0 ;
    for j = 1:n_time
       q_t = integrateOverdampedOscillators(q_t,k_dt_over_2_alpha,sigma_beta) ;
       q(:,j+1) = q_t ;
    end

end