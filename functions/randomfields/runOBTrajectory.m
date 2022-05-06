function f_t = runOBTrajectory(fluct_info,dt,n_steps,n_steps_fluct)
% computes the time-averaged Ornstein-Uhlenbeck process

% set up array for fluctuation averages
f_t = zeros([fluct_info.n_OB,n_steps]) ;
% get time step for fluctuations
dt_fluct = dt / n_steps_fluct ;
% sample initial fluctuation values
f_0 = sqrt(transpose(fluct_info.Deltaf_sq_OB)).*randn([fluct_info.n_OB,1]) ;
% set up arrays overdamped oscillator variables for integration
k_f = ones([fluct_info.n_OB,1])./transpose(fluct_info.Deltaf_sq_OB) ;
alpha = transpose(fluct_info.tau_OB).*k_f ;
% mean_f_t_sq = zeros(size(f_0)) ;
for k = 1:n_steps
    % run a segment of a trajectory
    f_t_segment = generateOverdampedOscillatorTrajectory(f_0, alpha, k_f, 1.0, dt_fluct, n_steps_fluct+1) ;
%     mean_f_t_sq = mean_f_t_sq + sum(f_t_segment.*f_t_segment,2) ;
    % average the segment
    f_t(:,k) = (dt_fluct/dt)*trapz(f_t_segment,2) ;
    % get new initial condition
    f_0 = f_t_segment(:,end) ;
end

% mean_f_t_sq = mean_f_t_sq / (n_steps*n_steps_fluct) ;
% mean_f_t_sq 

end