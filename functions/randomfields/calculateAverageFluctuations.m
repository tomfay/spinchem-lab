function omega_t = calculateAverageFluctuations(B_t_full,n_stoch,dt,n_fluct, n_time)

omega_t = zeros([n_fluct,n_time-1]) ;
t_stoch = (0:(n_stoch))*dt ;
for i = 1:(n_time-1)
    for m = 1:n_fluct
        omega_t(m,i) = (1.0/(n_stoch*dt))*trapz( t_stoch , B_t_full(m,(1+(i-1)*n_stoch):(i*n_stoch+1)) ) ;
    end
    %                     omega_t(1,i) = (1.0/(n_stoch*dt))*trapz( t_stoch , B_t_full(1,(1+(i-1)*n_stoch):(i*n_stoch+1)) ) ;
    %                     omega_t(2,i) = (1.0/(n_stoch*dt))*trapz( t_stoch , B_t_full(2,(1+(i-1)*n_stoch):(i*n_stoch+1)) ) ;
    %                     omega_t(3,i) = (1.0/(n_stoch*dt))*trapz( t_stoch , B_t_full(3,(1+(i-1)*n_stoch):(i*n_stoch+1)) ) ;
end


end