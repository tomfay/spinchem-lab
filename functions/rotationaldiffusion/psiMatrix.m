function psi = psiMatrix(q)

    
    psi = [ -0.5 * (q(2:4)' ) ;
                0.5 * ( q(1) * eye(3) - cross(repmat(q(2:4),[1,3]),eye(3)) ) ] ;


end