function deriv_q = derivativeQuaternion( q , omega )

    deriv_q = [ -0.5 * (q(2:4)' * omega) ;
                0.5 * ( q(1) * omega - cross( q(2:4) , omega ) ) ] ;

end