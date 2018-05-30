function [ r, s ] = assembleInpterpNode( Nh, Nz )
    [ xh, ~ ] = zwglj( Nh + 1 );
    [ xz, ~ ] = zwglj( Nz + 1 );
    r = xh * ones( 1, Nz + 1 );
    s = ones( Nh + 1, 1 ) * xz';
    
    r = r(:); 
    s = s(:);
end