function rotate_q = rotateQuaternionByOmega( q , omega, dt )

    omega_norm = norm(omega)  ;
    omega_normed = omega * (1.0/omega_norm) ;
    
    rotate_q = cos( omega_norm * dt * 0.5 ) * q ...
               + sin( omega_norm * dt * 0.5 ) * [ -  q(2:4)' * omega_normed ;
                 ( q(1) * omega_normed + cross( omega_normed , q(2:4) ) ) ] ;

end