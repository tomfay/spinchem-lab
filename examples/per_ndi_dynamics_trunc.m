% Per-NDI radical pair dynamics with 5 Per and 4 NDI HFCCs
% add the functions methods to the path (change as needed)
addpath(genpath('../functions'))
% constants
gamma_e = 176.0859644 ; % gyromagnetic ratio of free electron spin in mT^-1 mus^-1
g_e = 2.00231930436256 ; % free electron g-factor

% set k_S, k_T (in mus^-1)
k_S_permus = 20 ;
k_T_permus = 0.5 ;

% include triplet fraction, f_T = triplet fraction, set to true to include
include_triplet_fraction = true ;
f_T = 0.05 ;

% the important outputs are t, O_t_on, sigma_O_t_on, O_t_off, and
% sigma_O_t_off
% t is the time grid [1 x n_t]
% O_t is an array of <1>, <P_S>, <P_T+> , <P_T0>, <P_T-> at different times
% [5 x n_t]
% sigma_O_t is the standard error in each of the above 
% on/off refers to B = 25 mT or 0 mT

% g factors
g_per = (1/3)*(2.00218+2.00318+2.00238) ;
g_ndi = 2.0040 ;
% set up radical pair system
spin_system = struct() ;
spin_system.type = "radical pair" ;
spin_system.use_symmetry = true ;
spin_system.omega_1 = [0;0;0.0] ;
spin_system.omega_2 = [0;0;0.0] ;
spin_system.exchange = 0 ; % exchange term used is +2*exchange S_1.S_2
spin_system.dipolar_tensor = zeros([3,3]) ;
spin_system.A_1_tensors = kron([-0.410,-0.310],speye(3)) ;
spin_system.A_2_tensors = kron([-0.0963,-0.1927],speye(3)) ;
spin_system.N_1 = [3,2] ;
spin_system.N_2 = [2,2] ;
spin_system.g_1 = [2,2] ;
spin_system.g_2 = [3,2] ;
spin_system.k_S = k_S_permus/gamma_e ;
spin_system.k_T = k_T_permus/gamma_e ;

Z = prod(spin_system.g_1.^spin_system.N_1)*prod(spin_system.g_2.^spin_system.N_2) ;

% set up dynamics
dynamics = struct() ;
% sampling
dynamics.sampling = struct() ;
dynamics.sampling.method = "SU(N)" ;
dynamics.sampling.n_sample = 5 ;
dynamics.sampling.initial_electron_spin_state = "singlet" ;
% this option discards large Z blocks in the symmetry adapted calculation,
% the resulting observables need renormalisation in post processing
dynamics.sampling.symmetry_block_truncation = "discard large Z" ;
dynamics.sampling.truncation_tol = 1e-2 ;
% integration parameters - increase dt and n_steps to go to longer times
dynamics.integrator = struct() ;
dynamics.integrator.method = "adaptive SIA" ;
dynamics.integrator.dt = 2e-3*gamma_e ;
dynamics.integrator.n_steps = 1250 ;
dynamics.integrator.stride = 2 ;
dynamics.integrator.dim_krylov = 16 ;
dynamics.integrator.tol = 1e-7 ;
% observables
dynamics.observables = struct() ;
dynamics.observables.O_el = {speye(4),(singProj()),tripProjPlus(),(tripProj()-tripProjPlus()-tripProjMinus()),tripProjMinus()} ;

% run dynamics
tic
[O_t_off,t,sigma_O_t_off] = runQMSpinDynamics(spin_system,dynamics) ;
O_t_off = O_t_off * (1/O_t_off(1,1)) ;
sigma_O_t_off = sigma_O_t_off * (1/O_t_off(1,1)) ;
toc

% run dynamics with field on
spin_system.omega_1 = [0;0;25]*(g_per/g_e) ;
spin_system.omega_2 = [0;0;25]*(g_ndi/g_e) ;
dynamics.sampling.initial_electron_spin_state = "singlet" ;
tic
[O_t_on,t,sigma_O_t_on] = runQMSpinDynamics(spin_system,dynamics) ;
O_t_on = O_t_on * (1/O_t_on(1,1)) ;
sigma_O_t_on = sigma_O_t_on * (1/O_t_on(1,1)) ;
toc

if (include_triplet_fraction)
    dynamics.sampling.initial_electron_spin_state = "triplet" ;
    spin_system.omega_1 = [0;0;0]*(g_per/g_e) ;
    spin_system.omega_2 = [0;0;0]*(g_ndi/g_e) ;
    tic
    [O_t_T_off,t,sigma_O_t_T_off] = runQMSpinDynamics(spin_system,dynamics) ;
    O_t_T_off = O_t_T_off * (1/O_t_T_off(1,1)) ;
    sigma_O_t_T_off = sigma_O_t_T_off * (1/O_t_T_off(1,1)) ;
    toc
    spin_system.omega_1 = [0;0;25]*(g_per/g_e) ;
    spin_system.omega_2 = [0;0;25]*(g_ndi/g_e) ;
    tic
    [O_t_T_on,t,sigma_O_t_T_on] = runQMSpinDynamics(spin_system,dynamics) ;
    O_t_T_on = O_t_T_on * (1/O_t_T_on(1,1)) ;
    sigma_O_t_T_on = sigma_O_t_T_on * (1/O_t_T_on(1,1)) ;
    toc
    O_t_S_on = O_t_on ;
    sigma_O_t_S_on = sigma_O_t_on ;
    O_t_S_off = O_t_off ;
    sigma_O_t_S_off = sigma_O_t_off ;
    O_t_on = (1-f_T) * O_t_S_on + (f_T) * O_t_T_on ;
    O_t_off = (1-f_T) * O_t_S_off + (f_T) * O_t_T_off ;
    sigma_O_t_on = sqrt( (1-f_T)^2 * sigma_O_t_S_on.*sigma_O_t_S_on + (f_T)^2 * sigma_O_t_T_on.*sigma_O_t_T_on) ;
    sigma_O_t_off = sqrt( (1-f_T)^2 * sigma_O_t_S_off.*sigma_O_t_S_off + (f_T)^2 * sigma_O_t_T_off.*sigma_O_t_T_off) ;
end

t_mus = t /gamma_e ;
