function corr_A_t_B_0 = calculateCorrAB(A,B,n)

    % evaluates <A(t)B(0)> from a trajectory of A and B values
    corr_A_t_B_0 = zeros([1,n]) ;
    for k = 1:n
       corr_A_t_B_0(k) = mean(A(1:(end-k+1)).*B(k:end)) ;
    end
    

end