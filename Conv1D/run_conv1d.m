function run_conv1d

ne = [40, 60, 80, 100];
k = [1, 2, 3];
len = 2./ne;

Nmesh = numel(ne);
Ndeg = numel(k);
dofs = zeros(Nmesh, Ndeg);
for n = 1:Ndeg
    for m = 1:Nmesh
        dofs(m, n) = ne(m) * (k(n)+1);
    end
end

errInf = zeros(Nmesh, Ndeg);
err2 = zeros(Nmesh, Ndeg);
err1 = zeros(Nmesh, Ndeg);

for n = 1:Ndeg
    for m = 1:Nmesh
        conv = conv1d_advection( k(n), ne(m) );
        conv.RK45_solve;
        err2(m, n) = conv.norm_err2(conv.ftime);
        err1(m, n) = conv.norm_err1(conv.ftime);
        errInf(m, n) = conv.norm_errInf(conv.ftime);
    end
end

%marker = {'r-o', 'b-s', 'k-*', 'c-^'};
for n = 1:Ndeg
    subplot(1,3,1); plot(len, err1(:, n)); hold on;
    subplot(1,3,2); plot(len, err2(:, n)); hold on;
    subplot(1,3,3); plot(len, errInf(:, n)); hold on;
end

for n = 1:3
    subplot(1,3,n);
    box on; grid on;
    set(gca, 'XScale', 'log', 'YScale', 'log');
    lendstr = cell(Ndeg, 1);
    for m = 1:Ndeg
        lendstr{m} = ['p=', num2str(m)];
    end
    legend(lendstr, 'box', 'off');
end


for n = 1:Ndeg
    fprintf('\n==================deg = %d==================\n', n);
    t = convergence_table(len, err1(:, n), err2(:, n), errInf(:, n))
end

end

function t1 = convergence_table(len, err1, err2, errInf)
t1 = table;
if isrow(len)
    len = len';
end
t1.len = len;
t1.('err1') = err1(:);
t1.('a1') = get_ratio(len, err1);

t1.('err2') = err2(:);
t1.('a2') = get_ratio(len, err2);

t1.('errf') = errInf(:);
t1.('af') = get_ratio(len, errInf);
end

function a = get_ratio(len, err)
Nmesh = numel(len);

a = zeros(Nmesh, 1);
for m = 2:Nmesh
    scal_ratio = log2( len(m)/len(m-1) );
    a(m) = log2( err(m)/err(m-1) )./scal_ratio;
end
end
