% standard triangle cell
tol = 1e-10;

tri_3 = ndg_lib.std_cell.tri(3);
tri_4 = ndg_lib.std_cell.tri(4);
tri_5 = ndg_lib.std_cell.tri(5);

%% interpolation points
r3 = load('Coor_Test/r_3.cc');
s3 = load('Coor_Test/s_3.cc');

for i = 1:tri_3.Np
    assert( abs(tri_3.r(i) - r3(i)) <= tol)
    assert( abs(tri_3.s(i) - s3(i)) <= tol)
end

r4 = load('Coor_Test/r_4.cc');
s4 = load('Coor_Test/s_4.cc');

for i = 1:tri_4.Np
    assert( abs(tri_4.r(i) - r3(i)) <= tol)
    assert( abs(tri_4.s(i) - s3(i)) <= tol)
end

r5 = load('Coor_Test/r_5.cc');
s5 = load('Coor_Test/s_5.cc');

for i = 1:tri_5.Np
    assert( abs(tri_5.r(i) - r3(i)) <= tol)
    assert( abs(tri_5.s(i) - s3(i)) <= tol)
end

%% orthgonal function
% test of the orthgonal functions.
% the LGL quadrature formula is exact for polynomials of degree ¡Ü 2N ? 1.

for i = 1:tri_3.Np
    for j = 1:( tri_3.N*2 - 1 - i )
        temp = sum(tri_3.w.*tri_3.V(:,i).*tri_3.V(:,j));
        assert( abs(temp - (i==j))<tol )
    end
end

for i = 1:tri_4.Np
    for j = 1:( tri_4.N*2 - 1 - i )
        temp = sum(tri_4.w.*tri_4.V(:,i).*tri_4.V(:,j));
        assert( abs(temp - (i==j))<tol )
    end
end

for i = 1:tri_5.Np
    for j = 1:( tri_5.N*2 -1 - i )
        temp = sum(tri_5.w.*tri_5.V(:,i).*tri_5.V(:,j));
        assert( abs(temp - (i==j))<tol )
    end
end

