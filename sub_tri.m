N = 5;
tri = ndg_lib.std_cell.tri(N);

EToV = zeros(3, 0.5*(N+1)*N+0.5*(N-1)*N);
sk = 1;
for row = 1:N
    v1 = tri.Fmask(row, 3);
    v2 = tri.Fmask(row+1, 3);
    for kb = 1:row
        EToV(:, sk) = [v1, v2, v2+1]';
        v1 = v1+1;
        v2 = v2+1;
        sk = sk+1;
    end
    
    v1 = tri.Fmask(row, 3);
    v2 = tri.Fmask(row+1, 3);
    for kt = 1:row-1
        EToV(:, sk) = [v1, v2+1, v1+1]';
        v1 = v1+1;
        v2 = v2+1;
        sk = sk+1;
    end
end

patch('Vertices',[tri.r, tri.s], 'Faces',EToV',...
    'FaceColor', [0.8, 0.9, 1]);

