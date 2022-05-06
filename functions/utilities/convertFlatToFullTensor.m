function A = convertFlatToFullTensor(A_flat)
    A = zeros([3,3]) ;
    A(1,2) = A_flat(4) ;
    A(1,3) = A_flat(5) ;
    A(2,3) = A_flat(6) ;
    A = A + transpose(A) ;
    A(1,1) = A_flat(1) ;
    A(2,2) = A_flat(2) ;
    A(3,3) = A_flat(3) ;
end