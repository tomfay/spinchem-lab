function rotate_q = rotateQuaternionByOmegaSecondOrder( q , omega, dt )

    omega_normed = ( 1.0 / norm(omega) ) * omega ;
    
    rotate_q = (1.0 - 0.125*(omega'*omega)*dt*dt)*q + derivativeQuaternion( q , omega ) * dt  ;

end