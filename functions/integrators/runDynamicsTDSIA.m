function [O_t,t,psi_t] = runDynamicsTDSIA(psi_0,H_0,A,f_t,n_steps,dt,O,dim_krylov,stride)

% numbers of observables, time points and O_t array
n_obs = numel(O);
n_t = n_steps+1 ;
n_fluct = numel(A) ;

if nargin==9
else
    stride = 1 ;
end
n_data = 1 + floor(n_steps/stride) ;
% empty array for observables
O_t = zeros([n_obs,n_data]) ;
t = (0:(n_data-1)) * (dt*stride) ;
% set up initial state 
psi_t = psi_0 ;
% calculate initial observables
for j = 1:n_obs
    O_t(j,1) = psi_t' * (O{j}*psi_t) ;
end
data_index = 2 ;
for k = 1:n_steps
    
    % f_t is the average fluctuation over each time step
    f = f_t(:,k) ;
    
    % constructs the effective Hamiltonian
    %     H_eff = H_0 + kron( electronHamiltonian(omega_eff,omega_eff,0,zeros(3,3)),id_nuc) ;
    H_eff = H_0 ;
    for j = 1:n_fluct 
        H_eff = H_eff + f(j) * A{j} ;
    end
    
    % propagates the state using matrix exponentiation
    %     psi_t = propagateExpM(H_eff,psi_t, dt) ;
    
    % propagates the state using the Cayley approximant to the propagator
    %     psi_t = propagateCayley(H_eff,psi_t,dt) ;
    
    % propagates the state using the short iterative arnoldi integrator
    psi_t = propagateSIA(H_eff,psi_t, dt,dim_krylov) ;
    
        % compute observables in the krylov subspace
    if mod(k,stride)==0
        for j = 1:n_obs
            O_t(j,data_index) = psi_t' * (O{j}*psi_t) ;
        end
        data_index = data_index + 1 ;
    end
end

end