tol = 1e-10;

%% Line: N = 1
N = 1;
Lin = StdRegions.Line(N);
Dr_ext = load('Dr-Test/Dr_1.cc');

for i = 1:Lin.nNode
    assert( abs(Lin.Dr(i) - Dr_ext(i)) <= tol);
end


%% Line: N = 2
N = 2;
Lin = StdRegions.Line(N);
Dr_ext = load('Dr-Test/Dr_2.cc');

for i = 1:Lin.nNode
    assert( abs(Lin.Dr(i) - Dr_ext(i)) <= tol);
end

%% Line: N = 3
N = 3;
Lin = StdRegions.Line(N);
Dr_ext = load('Dr-Test/Dr_3.cc');

for i = 1:Lin.nNode
    assert( abs(Lin.Dr(i) - Dr_ext(i)) <= tol);
end
%% Line: N = 4
N = 4;
Lin = StdRegions.Line(N);
Dr_ext = load('Dr-Test/Dr_4.cc');

for i = 1:Lin.nNode
    assert( abs(Lin.Dr(i) - Dr_ext(i)) <= tol);
end
%% Line: N = 5
N = 5;
Lin = StdRegions.Line(N);
Dr_ext = load('Dr-Test/Dr_5.cc');

for i = 1:Lin.nNode
    assert( abs(Lin.Dr(i) - Dr_ext(i)) <= tol);
end