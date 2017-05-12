% test of the LGL node set with exact coordinate
tol = 10e-10;

%% Line: N=1
r_ext = load('Coor-test/r_1.cc');

% triangle shape
N = 1;
Lin = StdRegions.Line(N);
for i = 1:Lin.nNode
    assert( abs(Lin.r(i) - r_ext(i)) <= tol)
end

%% Line: N=2
r_ext = load('Coor-test/r_2.cc');

% triangle shape
N = 2;
Lin = StdRegions.Line(N);
for i = 1:Lin.nNode
assert( abs(Lin.r(i) - r_ext(i)) <= tol)
end

%% Line: N=3
r_ext = load('Coor-test/r_3.cc');
% triangle shape
N = 3;
Lin = StdRegions.Line(N);
for i = 1:Lin.nNode
assert( abs(Lin.r(i) - r_ext(i)) <= tol)
end
%% Line: N=4
r_ext = load('Coor-test/r_4.cc');
% triangle shape
N = 4;
Lin = StdRegions.Line(N);
for i = 1:Lin.nNode
assert( abs(Lin.r(i) - r_ext(i)) <= tol)
end
%% Line: N=5
r_ext = load('Coor-test/r_5.cc');
% triangle shape
N = 5;
Lin = StdRegions.Line(N);
for i = 1:Lin.nNode
assert( abs(Lin.r(i) - r_ext(i)) <= tol)
end