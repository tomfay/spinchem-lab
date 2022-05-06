function rank_2_components = calculateRank2SphericalTensorComponents(A) 

    % a 1 x 5 matrix of the rank 2 spherical tensor components of a tensor
    % A.
    rank_2_components = 1.0i*zeros(1,5) ;
    
    % m = -2 component
    rank_2_components(1) = 0.5 * (A(1,1) - A(2,2) + 1.0i * (A(1,2)+A(2,1)) ) ;
    
     % m = -1 component
    rank_2_components(2) = 0.5 * (A(1,3) + A(3,1) + 1.0i * (A(2,3)+A(3,2)) ) ;
    
     % m = 0 component
    rank_2_components(3) = (1.0/sqrt(6.0)) * (2.0 * A(3,3)- (A(1,1)+A(2,2)) ) ;
    
     % m = +1 component
    rank_2_components(4) = -0.5 * (A(1,3) + A(3,1) - 1.0i * (A(2,3)+A(3,2)) ) ;
    
    % m = +2 component
    rank_2_components(5) = 0.5 * (A(1,1) - A(2,2) - 1.0i * (A(1,2)+A(2,1)) ) ;

end