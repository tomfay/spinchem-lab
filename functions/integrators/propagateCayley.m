function psi_new = propagateCayley(H,psi_0,dt) 
    psi_new = (speye(size(H)) + (0.5i*dt)* H)\( (speye(size(H)) - (0.5i*dt)* H) * psi_0) ;
end