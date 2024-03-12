% Per-NDI radical pair dynamics with 5 Per and 4 NDI HFCCs
% add the functions methods to the path (change as needed)
addpath(genpath('../functions'))
% constants
gamma_e = 176.0859644 ; % gyromagnetic ratio of free electron spin in mT^-1 mus^-1
g_e = 2.00231930436256 ; % free electron g-factor

% set k_S, k_T (in mus^-1)
k_S_permus = 0 ;
k_T_permus = 0 ;

% include triplet fraction, f_T = triplet fraction, set to true to include
include_triplet_fraction = false ;
f_T = 0.00 ;

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
g_1 = 2.0031 ;
g_2 = 2.0034 ;
D_zz = 0.1265 ;
xi = pi/2 ;
D_tensor = (D_zz/2) * (1-3*cos(xi)^2) * diag([-1,-1,2]) ;
B_0 = 345 ; 
B_1 = 0.03 ;
Omega = 345 ;


% set up radical pair system
spin_system = struct() ;
spin_system.type = "radical pair" ;
spin_system.use_symmetry = true ;
spin_system.omega_1 = [B_1;0;(g_1/g_e)*B_0 - Omega ] ;
spin_system.omega_2 = [B_1;0;(g_2/g_e)*B_0 - Omega] ;
spin_system.exchange = 0.01 ; % exchange term used is +2*exchange S_1.S_2
spin_system.dipolar_tensor = D_tensor ;
spin_system.A_1_tensors = kron([0.1],spdiags([0;0;1],[0],3,3)) ;
spin_system.A_2_tensors = kron([0.0621],spdiags([0;0;1],[0],3,3)) ;
spin_system.N_1 = [4] ;
spin_system.N_2 = [2] ;
spin_system.g_1 = [2] ;
spin_system.g_2 = [3] ;
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
dynamics.integrator.n_steps = 200 ;
dynamics.integrator.stride = 2 ;
dynamics.integrator.dim_krylov = 16 ;
dynamics.integrator.tol = 1e-7 ;
% observables
dynamics.observables = struct() ;
dynamics.observables.O_el = {speye(4),(singProj()),tripProjPlus(),(tripProj()-tripProjPlus()-tripProjMinus()),tripProjMinus(),sumSpinOpsY([2,2],[1,1])} ;

% run dynamics
tic
DB_vals = -3:0.1:3 ;
n_fields = numel(DB_vals) ;
spec = zeros(ceil(dynamics.integrator.n_steps/dynamics.integrator.stride)+1,n_fields) ;
for i = 1:n_fields
rng('default')
spin_system.omega_1 = [B_1;0;(g_per/g_e)*(B_0+DB_vals(i)) - Omega] ;
spin_system.omega_2 = [B_1;0;(g_ndi/g_e)*(B_0+DB_vals(i)) - Omega] ;
[O_t,t,sigma_O_t] = runQMSpinDynamics(spin_system,dynamics) ;
O_t = O_t * (1/O_t(1,1)) ;
sigma_O_t = sigma_O_t * (1/O_t(1,1)) ;
spec(:,i) = real(O_t(end,:)) ;
toc
drawnow
end

DB_fine = -3:0.01:3 ;
n_fine = numel(DB_fine) ;
spec_fine = zeros(numel(t),n_fine) ;
n_window_B = 200 ;
for i = 1:numel(t)
spec_fine(i,:) = smoothdata(spline(DB_vals,spec(i,:),DB_fine),'gaussian',n_window_B) ;
end

spec_smooth = spec_fine ;
n_window_t = 10 ;
for i = 1:n_fine
spec_smooth(:,1) = smoothdata(spec_fine(:,i),'gaussian',n_window_t ) ;
end

t_mus = t /gamma_e ;

