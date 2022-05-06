function psi_new = propagateExpM(H,psi_0,dt) 
    psi_new = expm( (-1.0i*dt)* H ) * psi_0 ;
end