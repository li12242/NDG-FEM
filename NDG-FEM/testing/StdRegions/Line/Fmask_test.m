
tol = 10e-10;

%% Line: N=1
Fmask_ext = load('Fmask-test/Fmask_1.cc');

% triangle shape
N = 1;
Lin = StdRegions.Line(N);
for i = 1:numel(Fmask_ext)
    assert( abs(Lin.Fmask(i) - Fmask_ext(i)) <= tol)
end

%% Line: N=2
Fmask_ext = load('Fmask-test/Fmask_2.cc');

% triangle shape
N = 2;
Lin = StdRegions.Line(N);
for i = 1:numel(Fmask_ext)
assert( abs(Lin.Fmask(i) - Fmask_ext(i)) <= tol)
end

%% Line: N=3
Fmask_ext = load('Fmask-test/Fmask_3.cc');
% triangle shape
N = 3;
Lin = StdRegions.Line(N);
for i = 1:numel(Fmask_ext)
assert( abs(Lin.Fmask(i) - Fmask_ext(i)) <= tol)
end
%% Line: N=4
Fmask_ext = load('Fmask-test/Fmask_4.cc');
% triangle shape
N = 4;
Lin = StdRegions.Line(N);
for i = 1:numel(Fmask_ext)
assert( abs(Lin.Fmask(i) - Fmask_ext(i)) <= tol)
end
%% Line: N=5
Fmask_ext = load('Fmask-test/Fmask_5.cc');
% triangle shape
N = 5;
Lin = StdRegions.Line(N);
for i = 1:numel(Fmask_ext)
assert( abs(Lin.Fmask(i) -Fmask_ext(i)) <= tol)
end