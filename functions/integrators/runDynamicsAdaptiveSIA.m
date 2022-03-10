function [O_t, t, psi_t] = runDynamicsAdaptiveSIA(psi_0,H,n_steps,dt,O,dim_krylov,tol,stride)

%evolve the system evolving under
n_time = n_steps+1 ;
n_obs = length(O) ;
d = size(H,1) ;
if nargin==8
else
    stride = 1 ;
end
n_data = 1 + floor(n_steps/stride) ;
% empty array for observables
O_t = zeros([n_obs,n_data]) ;
t = (0:(n_data-1)) * (dt*stride) ;

% set up krylov subspace and full space states
c_t = zeros([dim_krylov,1]) ;
c_t(1) = norm(psi_0) ;
c_t_corr = zeros([dim_krylov+1,1]) ;
c_t_corr(1) = c_t(1) ;
psi_t = psi_0 ;

% generate initial krylov subspace
[H_krylov, krylov_basis] = generateKrylovSubspace(H,psi_t,dim_krylov+1) ;
% create the propagator in the krylov subspace
U_krylov_dt = expm((-1.0i * dt) * H_krylov(1:dim_krylov,1:dim_krylov)) ;
U_krylov_dt_corr = expm((-1.0i * dt) * H_krylov) ;
% set up empty observable operators in krylov subspace
O_krylov = cell([1,n_obs]) ;
for j = 1:n_obs
    O_krylov{j} = krylov_basis(:,1:dim_krylov)' * full(O{j} * krylov_basis(:,1:dim_krylov))  ;
end
% krylov subspace observable matrix
% O_krylov_mat = kron(speye(n_obs),krylov_basis') * O_mat * krylov_basis ;

% calculate initial observables
for j = 1:n_obs
    O_t(j,1) = psi_t' * (O{j}*psi_t) ;
end
data_index = 2 ;

for k = 1:n_steps
    % check to see if Krylov subspace needs to be re-generated
    if norm([c_t;0]-c_t_corr) > tol* norm(c_t)
        % generate initial krylov subspace
        psi_t = krylov_basis(:,1:dim_krylov) * c_t ;
        [H_krylov, krylov_basis] = generateKrylovSubspace(H,psi_t,dim_krylov+1) ;
        % create the propagator in the krylov subspace
        U_krylov_dt = expm((-1.0i * dt) * full(H_krylov(1:dim_krylov,1:dim_krylov))) ;
        U_krylov_dt_corr = expm((-1.0i * dt) * full(H_krylov)) ;
        % set up observable operators in krylov subspace
        for j = 1:n_obs
            O_krylov{j} = krylov_basis(:,1:dim_krylov)' * full(O{j} * krylov_basis(:,1:dim_krylov))  ;
        end
        %         O_krylov_mat = kron(speye(n_obs),krylov_basis') * O_mat * krylov_basis ;
        c_t = zeros([dim_krylov,1]) ;
        c_t(1) = norm(psi_t) ;
        c_t_corr = zeros([dim_krylov+1,1]) ;
        c_t_corr(1) = c_t(1) ;
    end

    % propagate the state in the krylov subspace
    c_t = U_krylov_dt * c_t ;
    c_t_corr = U_krylov_dt_corr * c_t_corr ;

    % compute observables in the krylov subspace
    if mod(k,stride)==0
        for j = 1:n_obs
            O_t(j,data_index) = c_t' * (O_krylov{j} * c_t) ;
        end
        data_index = data_index + 1 ;
    end
end

% compute final state vector
psi_t = krylov_basis(:,1:dim_krylov) * c_t ;

end