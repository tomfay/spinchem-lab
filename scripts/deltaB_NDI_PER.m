% calculate the anisotropic spherical tensor components of the HF tensors
A_NDI_flat{1} = [0.011586 0.032114 -0.043700 -0.101834 -0.000008 0.000014] ;
A_NDI_flat{2} = [0.011586 0.032114 -0.043700 0.101834 0.000014 0.000008];
A_NDI_flat{3} =[ 0.011586 0.032114 -0.043700 0.101834 0.000014 0.000008];
A_NDI_flat{4} =[ 0.011586 0.032114 -0.043700 -0.101834 -0.000008 0.000014];
A_NDI_flat{5} =[0.035200 0.034000 -0.069200 0.000000 0.000000 0.000010] ;
A_NDI_flat{6} = [ 0.035200 0.034000 -0.069200 0.000000 0.000000 0.000010];
I_NDI = [0.5,0.5,0.5,0.5,1,1] ;
B_hyp_aniso_NDI_sq = 0 ;
for k = 1:numel(A_NDI_flat)
    A_NDI{k} = convertFlatToFullTensor(A_NDI_flat{k}) ;
    A_NDI_spher{k} = calculateRank2SphericalTensorComponents(A_NDI{k})  ;
    B_hyp_aniso_NDI_sq = B_hyp_aniso_NDI_sq + sum(abs(A_NDI_spher{k}).^2)*I_NDI(k)*(I_NDI(k)+1) 
end
Delta_B_NDI = sqrt(B_hyp_aniso_NDI_sq/36) ;

A_Per{1} = diag(0.1*[1.62,-1.74,0.12]) ;
A_Per{2} = A_Per{1} ;
A_Per{3} = A_Per{1} ;
A_Per{4} = A_Per{1} ;
A_Per{5} = diag(0.1*[-0.24,0.26,-0.02]) ;
A_Per{6} = A_Per{5} ;
A_Per{7} = A_Per{5} ;
A_Per{8} = A_Per{5} ;
A_Per{9} = diag(0.1*[2.14,-2.29,0.15]) ;
A_Per{10} = A_Per{9} ;
A_Per{11} = A_Per{9} ;
I_Per = 0.5*ones([1,11]);
B_hyp_aniso_Per_sq = 0 ;
for k = 1:numel(A_Per)
    A_Per_spher{k} = calculateRank2SphericalTensorComponents(A_Per{k})  ;
    B_hyp_aniso_Per_sq = B_hyp_aniso_Per_sq + sum(abs(A_Per_spher{k}).^2)*I_Per(k)*(I_Per(k)+1) 
end
Delta_B_Per = sqrt(B_hyp_aniso_Per_sq/36) ;