function v_rotated = rotateQuaternion(q,v)

    v_rotated = q(2:4) * (q(2:4)' * v) + ...
        q(1) * cross(repmat(q(2:4),[1,size(v,2)]),v) + ...
        ( q(1) * q(1) - 0.5 ) * v ;
    v_rotated = 2.0 * v_rotated ;
        
end