function I_x = spinOperatorX(g_I)
    
%     switch g_I
%         case 2
%             I_x = [0.0,0.5; 0.5,0.0] ;
%         case 3
%             I_x = zeros(3) ;
%             I_x(1,2) = 1.0/sqrt(2.0) ;
%             I_x(2,1) = I_x(1,2) ;
%             I_x(2,3) = I_x(1,2) ;
%             I_x(3,2) = I_x(1,2) ;
%         otherwise
%             I_x = [0] ;
%     end
%     I_x = sparse(I_x) ;
    I = 0.5*(g_I-1) ;
    Ms = (I:(-1):(-I+1))' ;
    band_l = [0.5*sqrt(I*(I+1) - Ms.*(Ms-1));0] ;
    band_u = [0;0.5*sqrt(I*(I+1) - Ms.*(Ms-1))] ;
%     I_x = gallery('tridiag',band,zeros([1,g_I]),band) ;
    I_x = spdiags([band_l zeros([g_I,1]) band_u],-1:1,g_I,g_I);
end