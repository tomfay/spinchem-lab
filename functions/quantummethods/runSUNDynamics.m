function [O_t,sigma_O_t,t] = runSUNDynamics(H,Z,dynamics)  
    % dimensionality of system
    d = 4 * Z ;
    % psi_0s
    if (dynamics.sampling.initial_electron_spin_state == "singlet")
        psi_0s = [[0;1/sqrt(2);-1/sqrt(2);0]] ;
        psi_0_weights = [1] ;
    elseif (dynamics.sampling.initial_electron_spin_state == "triplet")
        psi_0s = [[1;0;0;0],[0;1/sqrt(2);1/sqrt(2);0],[0;0;0;1]] ;
        psi_0_weights = [1,1,1] ;
    end
    n_psi_0s = size(psi_0s,2) ;
    
    n_obs = numel(dynamics.observables.O) ;
    n_data = 1+floor(dynamics.integrator.n_steps/dynamics.integrator.stride) ;
    O_t = zeros([n_obs,n_data]) ;
    sigma_O_t = zeros([n_obs,n_data]) ;

    for r = 1:dynamics.sampling.n_sample
        psi_nuc_0 = createSUNCoherentState(Z) ;
        for n = 1:n_psi_0s 
            psi_0 = kron(psi_0s(:,n),psi_nuc_0) ;
            if (dynamics.integrator.method == "adaptive SIA")
                [O_t_traj, t] = runDynamicsAdaptiveSIA(psi_0,H,dynamics.integrator.n_steps,...
                    dynamics.integrator.dt,dynamics.observables.O,...
                    dynamics.integrator.dim_krylov,dynamics.integrator.tol,...
                    dynamics.integrator.stride) ;
            end
            O_t_sample = O_t_traj * psi_0_weights(n) ;
        end
        O_t = O_t + O_t_sample ;
        sigma_O_t = sigma_O_t + O_t_sample.*O_t_sample ;
    end

    O_t = (1/dynamics.sampling.n_sample)*O_t ;
    n_sample = dynamics.sampling.n_sample ;
    sigma_O_t = (1/sqrt(n_sample))*sqrt((n_sample/(n_sample-1))*((1/n_sample)*sigma_O_t - O_t.*O_t)) ;


end