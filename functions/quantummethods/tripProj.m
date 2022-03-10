function P_T = tripProj()

    s_z = sparse([[0.5,0];[0,-0.5]]) ;
    s_plus = sparse([[0,1];[0,0]]) ;
    s_minus = sparse([[0,0];[1,0]]) ;
    
    P_T = 0.75*speye(4) + kron( s_z,s_z) + kron(0.5*s_plus,s_minus) ...
        + kron(0.5*s_minus,s_plus) ;

end