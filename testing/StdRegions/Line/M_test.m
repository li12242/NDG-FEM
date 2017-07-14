tol = 1e-10;

%% Line: N = 1
N = 1;
Lin = StdRegions.Line(N);
M_ext = load('M-Test/M_1.cc');

for i = 1:Lin.nNode
    assert( abs(Lin.M(i) - M_ext(i)) <= tol);
end


%% Line: N = 2
N = 2;
Lin = StdRegions.Line(N);
M_ext = load('M-Test/M_2.cc');

for i = 1:Lin.nNode
    assert( abs(Lin.M(i) - M_ext(i)) <= tol);
end
%% Line: N = 3
N = 3;
Lin = StdRegions.Line(N);
M_ext = load('M-Test/M_3.cc');

for i = 1:Lin.nNode
    assert( abs(Lin.M(i) - M_ext(i)) <= tol);
end
%% Line: N = 4
N = 4;
Lin = StdRegions.Line(N);
M_ext = load('M-Test/M_4.cc');

for i = 1:Lin.nNode
    assert( abs(Lin.M(i) - M_ext(i)) <= tol);
end
%% Line: N = 5
N = 5;
Lin = StdRegions.Line(N);
M_ext = load('M-Test/M_5.cc');

for i = 1:Lin.nNode
    assert( abs(Lin.M(i) - M_ext(i)) <= tol);
end

