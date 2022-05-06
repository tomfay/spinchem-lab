function [q_new] = integrateOverdampedOscillators(q_old,k_dt_over_2_alpha,sigma_beta)
    % q_old - the old coordinates
    % k_dt_over_alpha - array of k_n * dt / (2 alpha_n), where these are the
    % force constants, time step and friction coefficient of each
    % oscillator.
    % sigma_beta - array of variances of the random forces on each
    % oscillator, sigma_n = sqrt(2 *k_B * T * dt/alpha_n)


    % generate the random forces integrated over dt
    % variance is 2 alpha k_B T dt
    beta_new = sigma_beta .* randn(size(q_old)) ;
    
    % generate the new q values from the old and random force
    % uses trapezium rule for (1/alpha)integral F(q) dt and uses explicit
    % solution.
    q_new = ( (1-k_dt_over_2_alpha).*q_old + beta_new)./(1+k_dt_over_2_alpha) ;

end