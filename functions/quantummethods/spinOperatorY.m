function I_y = spinOperatorY(g_I)
    
%     switch g_I
%         case 2
%             I_y = [0.0,-0.5i; 0.5i,0.0] ;
%         case 3
%             I_y = 1.0i*zeros(3) ;
%             I_y(1,2) = -1.0i/sqrt(2.0) ;
%             I_y(2,1) = -I_y(1,2) ;
%             I_y(2,3) = I_y(1,2) ;
%             I_y(3,2) = -I_y(1,2) ;
%         otherwise 
%             I_y = [0] ;
%     end
%     I_y = sparse(I_y) ;
    I = 0.5*(g_I-1) ;
    Ms = (I:(-1):(-I+1))' ;
    band_l = [0.5i*sqrt(I*(I+1) - Ms.*(Ms-1));0] ;
    band_u = [0; -0.5i*sqrt(I*(I+1) - Ms.*(Ms-1))] ;
%     I_y = gallery('tridiag',band,zeros([1,g_I]),-band) ;
    I_y = spdiags([band_l zeros([g_I,1]) band_u],-1:1,g_I,g_I);
end