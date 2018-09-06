function [ mesh ] = makeUniformMesh1d( N, xlim, K, bcType )

% Nmesh = checkMeshNum( N, xtick, Mx );
% [ Nv, vx, EToV, EToR ] = makeMesh( Nmesh, Mx );

% meshLine = cell(Nmesh, 1);
% for m = 1:Nmesh
%     stdcell = StdLine( N(m) );
%     meshLine{m} = NdgMesh1d(stdcell, Nv, vx, Mx(m), EToV{m}, EToR{m});
% end
Nv = K + 1;
vx = linspace( xlim(1), xlim(2), Nv )';
EToV = [ 1:(Nv - 1); 2:Nv ]; 
EToR = enumRegion.Normal * ones( K, 1, 'int8' );

stdcell = StdLine( N );
mesh = NdgMesh1d(stdcell, Nv, vx, K, EToV, EToR);
% mesh = makeMeshUnion( Nmesh, meshLine{:} );
end

% function [ Nv, vx, EToV, EToR ] = makeMesh( Nmesh, Mx, xlim )
% Nv = sum(Mx) + 1;
% vx = linspace( xlim(1), xlim(2), Mx(1)+1 )';
% % for m = 2:Nmesh
% %     x = linspace( xlim(m), xlim(m+1), Mx(m)+1 );
% %     vx = [ vx; x(2:end)' ];
% % end
% % BCToV = [1, Nv; double( bcType(1) ), double( bcType(2) ) ];
% 
% EToV = cell( Nmesh, 1 );
% EToR = cell( Nmesh, 1 );
% for m = 1:Nmesh
%     K = Mx(m);
% 
%     finalVertId = sum( Mx(1:m) ) + 1;
%     startVertId = finalVertId - K;
%     EToV{m} = [ startVertId:(finalVertId-1); 
%         (startVertId+1):finalVertId ];
%     EToR{m} = enumRegion.Normal * ones( K, 1, 'int8' );
% end
% 
% end
% 
% function Nmesh = checkMeshNum( N, xtick, Mx )
% 
% Nmesh = numel( N );
% 
% if Nmesh ~= numel( Mx )
%     msgID = [ mfilename, ':InputMeshNumberError'];
%     msgtext = 'The number of input Mx is not equal to N.';
%     throw( MException(msgID, msgtext) );
% end
% 
% if Nmesh ~= ( numel( xtick ) - 1 )
%     msgID = [ mfilename, ':InputLengthNumberError'];
%     msgtext = 'The number of input mesh region xtick is not equal to N.';
%     throw( MException(msgID, msgtext) );
% end
% 
% % check the xtick is monotony
% dx = diff( xtick );
% if all( dx > 0 ) || all( dx < 0 )
%     return;
% else
%     msgID = [ mfilename, ':InputXTickError'];
%     msgtext = 'The input vector xtick should be monotony.';
%     throw( MException(msgID, msgtext) );
% end
% 
% end
