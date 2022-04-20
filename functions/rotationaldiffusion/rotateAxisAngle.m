function v_rotated = rotateAxisAngle(axis, angle, v)

    cos_angle = cos(angle) ;
    v_rotated = cos_angle * v + sin(angle) * cross(repmat(axis,[1,size(v,2)]),v) + axis * (1.0 - cos_angle) * (axis' * v)  ;

end