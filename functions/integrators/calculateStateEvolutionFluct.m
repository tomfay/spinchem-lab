function O_t = calculateStateEvolutionFluct(H_0,f_t,A,t,psi_0,O)

% numbers of observables, time points and O_t array
n_O = numel(O);
n_time = length(t) ;
O_t = (1.0+0.0i)*zeros(n_O,n_time) ;
n_fluct = numel(A) ;

% set up initial state 
psi_t = psi_0 ;
for k = 1:(n_time-1)
    % calculate the observables at the current time
    for i = 1:n_O
        O_t(i,k) = psi_t' * ( O{i} * psi_t ) ;
    end
    
    % calculate the average omega over the time step from t(k) to t(k+1)
    %     omega_eff = 0.5*(omega_t(:,k) + omega_t(:,k+1)) ;
    f = f_t(:,k) ;
    
    % constructs the effective Hamiltonian
    %     H_eff = H_0 + kron( electronHamiltonian(omega_eff,omega_eff,0,zeros(3,3)),id_nuc) ;
    H_eff = H_0 ;
    for j = 1:n_fluct 
        H_eff = H_eff + f(j) * A{j} ;
    end
    
    % calculate the time interval (not necessarily constant)
    dt = t(k+1) - t(k) ;
    
    % propagates the state using matrix exponentiation
    %     psi_t = propagateExpM(H_eff,psi_t, dt) ;
    
    % propagates the state using the Cayley approximant to the propagator
    %     psi_t = propagateCayley(H_eff,psi_t,dt) ;
    
    % propagates the state using the short iterative arnoldi integrator
    psi_t = propagateSIA(H_eff,psi_t, dt,4) ;
    
end

% calculate the final observables
k = n_time ;
for i = 1:n_O
    O_t(i,k) = psi_t' * ( O{i} * psi_t ) ;
end

end