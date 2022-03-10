% test for dynamics
% constants
gamma_e = 176.0859644 ; % gyromagnetic ratio of free electron spin in mT^-1 mus^-1
g_e = 2.00231930436256 ; % free electron g-factor

% g factors
g_per = (1/3)*(2.00218+2.00318+2.00238) ;
g_ndi = 2.0040 ;

% set up radical pair system
spin_system = struct() ;
spin_system.type = "radical pair" ;
spin_system.use_symmetry = true ; 
spin_system.omega_1 = [0;0;0.0] ;
spin_system.omega_2 = [0;0;0.0] ;
spin_system.exchange = -0.25*0 ;
spin_system.dipolar_tensor = zeros([3,3]) ;
spin_system.A_1_tensors = kron([-0.410],speye(3)) ;
% spin_system.A_1_tensors = kron([-0.410,-0.410,-0.410],speye(3)) ;
spin_system.A_2_tensors = kron([-0.0963],speye(3)) ;
spin_system.N_1 = [3] ;
spin_system.N_2 = [2] ;
spin_system.g_1 = [2] ;
spin_system.g_2 = [3] ;
spin_system.k_S = 20/gamma_e ;
spin_system.k_T = 0.5/gamma_e ;
% spin_system.k_S = 0.5/gamma_e ;
% spin_system.k_T = 20/gamma_e ;
Z = prod(spin_system.g_1)*prod(spin_system.g_2) ;

% set up dynamics
dynamics = struct() ;
% sampling
dynamics.sampling = struct() ;
dynamics.sampling.method = "SU(N)" ;
dynamics.sampling.n_sample = 5 ;
dynamics.sampling.initial_electron_spin_state = "singlet" ;
% integration parameters
dynamics.integrator = struct() ;
dynamics.integrator.method = "adaptive SIA" ;
dynamics.integrator.dt = 2e-3*gamma_e ;
dynamics.integrator.n_steps = 1250 ;
dynamics.integrator.stride = 2 ;
dynamics.integrator.dim_krylov = 16 ;
dynamics.integrator.tol = 1e-7 ;
% observables
dynamics.observables = struct() ;
dynamics.observables.O = {kron(singProj(),speye(Z)),kron(tripProj(),speye(Z)),speye(4*Z)} ;

% run dynamics
% tic
% [O_t_off,t,sigma_O_t_off] = runQMSpinDynamics(spin_system,dynamics) ;
% toc
dynamics.sampling.initial_electron_spin_state = "triplet" ;
tic
[O_t_T_off,t,sigma_O_t_T_off] = runQMSpinDynamics(spin_system,dynamics) ;
toc

% run dynamics with field on
spin_system.omega_1 = [0;0;25]*(g_per/g_e) ;
spin_system.omega_2 = [0;0;25]*(g_ndi/g_e) ;
dynamics.sampling.initial_electron_spin_state = "singlet" ;
% tic
% [O_t_on,t,sigma_O_t_on] = runQMSpinDynamics(spin_system,dynamics) ;
% toc
% dynamics.sampling.initial_electron_spin_state = "triplet" ;
% tic
% [O_t_T_on,t,sigma_O_t_T_on] = runQMSpinDynamics(spin_system,dynamics) ;
% toc
