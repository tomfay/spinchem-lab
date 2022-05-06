function [O_t,t,sigma_O_t] = runQMSpinDynamics(spin_system,dynamics)

if (spin_system.type == "radical pair" & ~isfield(spin_system,"use_symmetry"))
    if (isfield(spin_system,'include_fluctuations'))
        if (spin_system.include_fluctuations)
            error("Including fluctuations is not currently implemented for this type of simulation.") ;
        end
    end
    % get dimensionality of spin system
    Z = prod(spin_system.g_1) * prod(spin_system.g_2) ;
    d = 4 * Z ;
    % construct the spin Hamiltonian
    H = electronHamiltonianQMSp(spin_system.omega_1 , spin_system.omega_2, ...
        spin_system.exchange, spin_system.dipolar_tensor,  ...
        spin_system.A_1_tensors, spin_system.A_2_tensors, ...
        spin_system.g_1,spin_system.g_2) ;
    % construct the reaction operator
    k_S = spin_system.k_S ;
    k_T = spin_system.k_T ;
    K = kron((0.5*k_S)*singProj()+(0.5*k_T)*tripProj(),speye(Z))  ;
    % construct the effective Hamiltonian
    H_eff = H - 1.0i*K ;
    if (dynamics.sampling.method == "SU(N)")
        [O_t,sigma_O_t,t] = runSUNDynamics(H_eff,Z,dynamics)  ;
    end
    if (dynamics.sampling.method == "full trace")
        [O_t,t] = runFullTraceDynamics(H_eff,Z,dynamics)  ;
        sigma_O_t = zeros(size(O_t)) ;
    end
    

elseif (spin_system.type == "radical pair" & spin_system.use_symmetry == true)
    
    % check to see if fluctuations are included & what type are included
    if (isfield(spin_system,'include_fluctuations'))
        if (spin_system.include_fluctuations)
            include_fluctuations = true ;
            if ( isfield(spin_system,'fluct_operators_OB_el') )
            include_fluct_OB_el = true ;
            include_fluct_OB = true ;
            n_fluct_el = numel(spin_system.fluct_operators_OB_el) ;
            end
        else
            include_fluctuations = false ;
        end
    else
        include_fluctuations = false ;
    end


    % if using symmetry first construct the symmetry blocks
    n_blocks_1 = numel(spin_system.N_1) ;
    n_blocks_2 = numel(spin_system.N_2) ;
    n_blocks = n_blocks_1 + n_blocks_2 ;
    % generate the subspaces for the different symmetry blocks
    [block_weights_1,g_blocks_1,n_subspaces_1] = createSymmetryBlocks(spin_system.g_1,spin_system.N_1) ;
    [block_weights_2,g_blocks_2,n_subspaces_2] = createSymmetryBlocks(spin_system.g_2,spin_system.N_2) ;
    block_weights = [block_weights_1,block_weights_2] ;
    g_blocks = [g_blocks_1,g_blocks_2] ;
    n_subspaces = [n_subspaces_1,n_subspaces_2]  ;
    % Find the total weights of each subspace
    [total_weights,total_weight_Zs, Z_subspaces] = generateSubspaceTotalWeights(block_weights,n_subspaces,g_blocks) ;

    Z = prod([spin_system.g_1.^spin_system.N_1,spin_system.g_2.^spin_system.N_2]) ;
    N_subspaces_tot = prod(n_subspaces) ;
    % determine subspaces to exclude 
    % if discarding large Z - discard the largest Z subspaces up to some
    % tolerance
    if (dynamics.sampling.symmetry_block_truncation == "discard large Z")
        [Z_subspaces_sort,sort_indices] = sort(Z_subspaces,'descend') ;
        l = 1 ;
        discarded_weight = total_weight_Zs(sort_indices(l));
        discarded_subspace_indices = [] ;
        while (discarded_weight < dynamics.sampling.truncation_tol*Z)
            Z_subspaces_sort(l) ;
            discarded_subspace_indices = [discarded_subspace_indices,sort_indices(l)] ;
            l = l + 1 ;
            discarded_weight = discarded_weight + total_weight_Zs(sort_indices(l)) ;
        end
    else
        discarded_subspace_indices = [] ;
    end

    % Run the dynamcis for each subspace and weight it appropriately
    n_data = 1+floor(dynamics.integrator.n_steps/dynamics.integrator.stride) ;
    n_obs = numel(dynamics.observables.O_el) ;
    O_t = zeros([n_obs,n_data]) ;
    sigma_O_t = zeros(size(O_t)) ;
    for J = 1:N_subspaces_tot
        if (~ismember(J,discarded_subspace_indices))
            % get multiindex
            j = getMultiIndexFromIndex(J,n_subspaces) ;
            g_1_J = zeros([1,n_blocks_1]) ;
            g_2_J = zeros([1,n_blocks_2]) ;
            for k = 1:n_blocks_1
                g_1_J(k) = g_blocks{k}(j(k)) ;
            end
            for k = (1):(n_blocks_2)
                g_2_J(k) = g_blocks{n_blocks_1+k}(j(n_blocks_1+k)) ;
            end
            H = electronHamiltonianQMSp(spin_system.omega_1 , spin_system.omega_2, ...
                spin_system.exchange, spin_system.dipolar_tensor,  ...
                spin_system.A_1_tensors, spin_system.A_2_tensors, ...
                g_1_J,g_2_J) ;
            K = kron((0.5*spin_system.k_S)*singProj()+(0.5*spin_system.k_T)*tripProj(),speye(Z_subspaces(J)))  ;
            H_eff = H - 1.0i*K ;
            dynamics_J = dynamics ;
            dynamics_J.O = {} ;
            for k = 1:n_obs
                dynamics_J.observables.O{k} = kron(dynamics.observables.O_el{k},speye(Z_subspaces(J))) ;
            end
            if (~include_fluctuations)
                % calculate the exact trace if the number of SU(N) samples
                % is greater than the Z of the subspace, otherwise use
                % SU(N) sampling
                if (Z_subspaces(J) > dynamics.sampling.n_sample)
                    [O_t_subspace,sigma_O_t_subspace,t] = runSUNDynamics(H_eff,Z_subspaces(J),dynamics_J)  ;
                    O_t = O_t + total_weight_Zs(J) * O_t_subspace ;
                    sigma_O_t = sigma_O_t + (total_weight_Zs(J)*total_weight_Zs(J)) * (sigma_O_t_subspace.*sigma_O_t_subspace) ;
                else
                    [O_t_subspace,t] = runFullTraceDynamics(H_eff,Z_subspaces(J),dynamics_J)  ;
                    O_t = O_t + total_weight_Zs(J) * O_t_subspace ;
                end
            else
                % create the set of fluctuation operators & information
                % about fluctuating terms int he Hamiltonian
                fluct_info = struct() ;
                fluct_info.A = {} ;
                if (include_fluct_OB_el)
                    for k = 1:n_fluct_el
                        fluct_info.A = [fluct_info.A,{kron(spin_system.fluct_operators_OB_el{k},speye(Z_subspaces(J)))}] ;
                    end
                    fluct_info.Deltaf_sq_OB = spin_system.fluct_Deltaf_sq_OB_el ;
                    fluct_info.tau_OB = spin_system.fluct_tau_OB_el ;
                    fluct_info.n_OB = numel(fluct_info.A) ;
                end
                % run the dynamics with fluctuations
                [O_t_subspace,sigma_O_t_subspace,t] = runSUNDynamicsTDFluct(H_eff,Z_subspaces(J),dynamics_J,fluct_info)  ;
                O_t = O_t + total_weight_Zs(J) * O_t_subspace ;
                sigma_O_t = sigma_O_t + (total_weight_Zs(J)*total_weight_Zs(J)) * (sigma_O_t_subspace.*sigma_O_t_subspace) ;
            end
            
        end
    end
    O_t = (1/Z)*O_t ;
    sigma_O_t = sqrt((1/(Z*Z))*sigma_O_t) ;

end


end