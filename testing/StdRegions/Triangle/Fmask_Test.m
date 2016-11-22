tol = 1e-10;

%% Tri: N = 3
N = 3;
tri = StdRegions.Triangle(N);
facelist = [1,2,3,4; 4, 7, 9, 10; 10, 8, 5, 1]';
for i = 1:tri.nFaceNode
    assert( abs(tri.Fmask(i) - facelist(i)) <= tol);
end% for

%% Tri: N = 4
N = 4;
tri = StdRegions.Triangle(N);
facelist = [1,2,3,4,5; 5,9,12,14,15;15,13,10,6,1]';
for i = 1:tri.nFaceNode
    assert( abs(tri.Fmask(i) - facelist(i)) <= tol);
end% for
%% Tri: N = 5
N = 5;
tri = StdRegions.Triangle(N);
facelist = [1,2,3,4,5,6; 6,11,15,18,20,21; 21,19,16,12,7,1]';
for i = 1:tri.nFaceNode
    assert( abs(tri.Fmask(i) - facelist(i)) <= tol);
end% for