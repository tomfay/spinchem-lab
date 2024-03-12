% Per-NDI radical pair dynamics with all isotropic HFCCs
% add the functions methods to the path (change as needed)
addpath(genpath('../functions'))
% constants
gamma_e = 176.0859644 ; % gyromagnetic ratio of free electron spin in mT^-1 mus^-1
g_e = 2.00231930436256 ; % free electron g-factor

% set k_S, k_T (in mus^-1)
k_S_permus = 0 ;
k_T_permus = 0 ;


% the important outputs are t, O_t_on, sigma_O_t_on, O_t_off, and
% sigma_O_t_off
% t is the time grid [1 x n_t]
% O_t is an array of <1>, <P_S>, <P_T+> , <P_T0>, <P_T-> at different times
% [5 x n_t]
% sigma_O_t is the standard error in each of the above quantities
% on/off refers to B = 50 muT at 0 or 90 degrees

% g factors
g_fad = g_e ;
g_trp = g_e ;
% set up radical pair system
spin_system = struct() ;
spin_system.type = "radical pair" ;
spin_system.use_symmetry = true ;
spin_system.omega_1 = [0;0;0.050]*(g_fad/g_e) ;
spin_system.omega_2 = [0;0;0.050]*(g_trp/g_e) ;
spin_system.exchange = 0 ; % exchange term used is +2*exchange S_1.S_2
spin_system.dipolar_tensor = zeros([3,3]) ;
A_N5 = [-0.0989,0.0039,0;
        0.0039,-0.0881,0;
        0,0,1.7569] ;
A_N10 = [-0.0190,-0.0048,0;
        -0.0048,-0.0196,0;
        0,0,0.6046] ;
spin_system.A_1_tensors = [A_N5,A_N10] ;
spin_system.A_2_tensors = [] ;
spin_system.N_1 = [1,1] ;
spin_system.N_2 = [] ;
spin_system.g_1 = [3,3] ;
spin_system.g_2 = [] ;
spin_system.k_S = k_S_permus/gamma_e ;
spin_system.k_T = k_T_permus/gamma_e ;

Z = prod(spin_system.g_1.^spin_system.N_1)*prod(spin_system.g_2.^spin_system.N_2) ;

% set up dynamics
dynamics = struct() ;
% sampling
dynamics.sampling = struct() ;
dynamics.sampling.method = "SU(N)" ;
dynamics.sampling.n_sample = 9 ; % increase this to get better statistics
dynamics.sampling.initial_electron_spin_state = "singlet" ;
% this option discards large Z blocks in the symmetry adapted calculation,
% the resulting observables need renormalisation in post processing
dynamics.sampling.symmetry_block_truncation = "discard large Z" ;
dynamics.sampling.truncation_tol = 1e-2 ;
% integration parameters - increase dt and n_steps to go to longer times
dynamics.integrator = struct() ;
dynamics.integrator.method = "adaptive SIA" ;
dynamics.integrator.dt = 0.1e-3*gamma_e ;
dynamics.integrator.n_steps = 10000 ;
dynamics.integrator.stride = 2 ;
dynamics.integrator.dim_krylov = 16 ;
dynamics.integrator.tol = 1e-10 ;
% observables
dynamics.observables = struct() ;
dynamics.observables.O_el = {speye(4),(singProj()),tripProjPlus(),(tripProj()-tripProjPlus()-tripProjMinus()),tripProjMinus()} ;

% run dynamics
tic
[O_t_off,t,sigma_O_t_off,info_off] = runQMSpinDynamics(spin_system,dynamics) ;
O_t_off = O_t_off * (1/O_t_off(1,1)) ;
sigma_O_t_off = sigma_O_t_off * (1/O_t_off(1,1)) ;
toc

% run dynamics with field on
spin_system.omega_1 = [0.05;0.0;0]*(g_fad/g_e) ;
spin_system.omega_2 = [0.05;0.0;0]*(g_trp/g_e) ;
dynamics.sampling.initial_electron_spin_state = "singlet" ;
tic
[O_t_on,t,sigma_O_t_on,info_on] = runQMSpinDynamics(spin_system,dynamics) ;
O_t_on = O_t_on * (1/O_t_on(1,1)) ;
sigma_O_t_on = sigma_O_t_on * (1/O_t_on(1,1)) ;
toc

% dynamics_full = dynamics ;
% dynamics_full.oberseables.O = {} ;
% for r = 1:numel(dynamics.observables.O_el)
%     dynamics_full.observables.O{r} = kron(dynamics.observables.O_el{r},speye(Z)) ;
% end
% [O_t,t] = runFullTraceDynamics(info_off.H,Z,dynamics_full)  ;

rho_0 = full(kron(singProj(),speye(Z))) ;
P_S = rho_0 ;
rho_0 = rho_0 / Z ;
U = expm(-1.0i * info_off.H*t(2)-t(1)) ;
U_dag = U' ;
O_t_full = zeros([1,numel(t)]) ;
rho_t = rho_0 ;
O_t_full(1) = sum(P_S.*rho_t,'all') ;
for r = 2:numel(t) 
    rho_t = U * rho_t * U_dag ;
    O_t_full(r) = sum(P_S.*rho_t,'all') ;
end

t_mus = t/gamma_e ;

