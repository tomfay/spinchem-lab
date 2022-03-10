function I_z = spinOperatorZ(g_I)
    
%     switch g_I
%         case 2
%             I_z = [0.5,0.0; 0.0,-0.5] ;
%         case 3
%             I_z = zeros(3) ;
%             I_z(1,1) = 1.0 ;
%             I_z(3,3) = -1.0 ;
%         otherwise
%             I_z = [0] ;
%     end
%     I_z = sparse(I_z) ;
    I = 0.5*(g_I-1) ;
    Ms = [I:(-1):(-I)]' ;
    I_z = spdiags([Ms],0,g_I,g_I) ;
end