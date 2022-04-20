function q = quaternionFromAxisAngle( axis , angle )

    q = [ cos(angle*0.5) ;
          sin(angle*0.5) * axis ] ;

end