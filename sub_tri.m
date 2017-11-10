N = 2;
% tri = ndg_lib.std_cell.tri(N);
tri = StdSubTri(N);
% EToV = zeros(3, 0.5*(N+1)*N+0.5*(N-1)*N);
% sk = 1;
% for row = 1:N
%     v1 = tri.Fmask(row, 3);
%     v2 = tri.Fmask(row+1, 3);
%     for kb = 1:row
%         EToV(:, sk) = [v1, v2, v2+1]';
%         v1 = v1+1;
%         v2 = v2+1;
%         sk = sk+1;
%     end
%     
%     v1 = tri.Fmask(row, 3);
%     v2 = tri.Fmask(row+1, 3);
%     for kt = 1:row-1
%         EToV(:, sk) = [v1, v2+1, v1+1]';
%         v1 = v1+1;
%         v2 = v2+1;
%         sk = sk+1;
%     end
% end

patch('Vertices',[tri.r, tri.s], 'Faces',tri.EToV',...
    'FaceColor', [0.8, 0.9, 1]);
hold on;

Nface = tri.Nface;
xp = zeros( tri.Nface, tri.NFV );
yp = zeros( tri.Nface, tri.NFV );
for n = 1:tri.NFV
    nodeId = tri.EToV(:, n);
    xc = mean( tri.r( nodeId ) );
    yc = mean( tri.s( nodeId ) );
    
    for f1 = 1:Nface
        f2 = mod(f1, Nface)+1;
        fv1 = tri.EToV(f1,n);
        fv2 = tri.EToV(f2,n);
        
        xf = ( tri.r(fv1) + tri.r(fv2) ) * 0.5;
        yf = ( tri.s(fv1) + tri.s(fv2) ) * 0.5;
        
        xp(f1, n) = xf;
        yp(f1, n) = yf;
        plot( [xc, xf], [yc, yf], 'r--' );
    end
end

quiver(xp(:), yp(:), tri.nr, tri.ns);
axis equal