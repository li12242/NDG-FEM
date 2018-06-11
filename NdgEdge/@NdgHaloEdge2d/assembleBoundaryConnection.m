function [ obj ] = assembleBoundaryConnection( obj, BCToV )
%ASSEMBLEBOUNDARYCONNECTION Summary of this function goes here
%   Detailed explanation goes here

Ne = obj.Ne;
Nb = size( BCToV, 2 );

ftype = zeros( Ne, 1 );
for n = 1:Ne
    vert = sort( obj.FToV(:, n) );
    
    for m = 1:Nb
        vert1 = sort( BCToV(1:2, m) );
        
        if all( vert == vert1 )
            ftype(n) = BCToV( 3, m );
            break;
        end
    end
end

obj.ftype = enumBoundaryCondition( ftype );

end

