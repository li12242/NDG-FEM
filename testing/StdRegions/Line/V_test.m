% test of Vandermonde matrix
tol = 1e-10;
%% Line N=1
V_ext = load('V-Test/V_1.cc');

% triangle shape
N = 1;
Lin= StdRegions.Line(N);
for i = 1:numel(V_ext)
    assert( abs(Lin.VandMatrix(i) - V_ext(i)) <= tol);
end

%% Line N=2
V_ext = load('V-Test/V_2.cc');

% triangle shape
N = 2;
Lin= StdRegions.Line(N);
for i = 1:numel(V_ext)
    assert( abs(Lin.VandMatrix(i) - V_ext(i)) <= tol);
end
%% Line N=3
V_ext = load('V-Test/V_3.cc');
% triangle shape
N = 3;
Lin= StdRegions.Line(N);

for i = 1:numel(V_ext)
    assert( abs(Lin.VandMatrix(i) - V_ext(i)) <= tol);
end
%% Line N=4
V_ext = load('V-Test/V_4.cc');

% triangle shape
N = 4;
Lin= StdRegions.Line(N);
for i = 1:numel(V_ext)
    assert( abs(Lin.VandMatrix(i) - V_ext(i)) <= tol);
end
%% Line N=5
V_ext = load('V-Test/V_5.cc');

% triangle shape
N = 5;
Lin= StdRegions.Line(N);
for i = 1:numel(V_ext)
    assert( abs(Lin.VandMatrix(i) - V_ext(i)) <= tol);
end