function run_pb2d

order = [1];
ne = [20, 40, 60];
len = 1./ne;
type = NdgCellType.Tri;

Nmesh = numel(ne);
Ndeg = numel(order);
dofs = zeros(Nmesh, Ndeg);
time = zeros(Nmesh, Ndeg);
for n = 1:Ndeg
    for m = 1:Nmesh
        dofs(m, n) = ne(m) * (order(n)+1).^2;
    end
end

errInf = zeros(Nmesh, Ndeg);
err2 = zeros(Nmesh, Ndeg);
err1 = zeros(Nmesh, Ndeg);

linewidth = 1.5; 
markersize = 8;
% c = 'b';
c = 'r';
% c = 'g';
% c = 'c';
% linestyle = '-';
linestyle = ':';
% linestyle = '-.';
% linestyle = '--';

color = {c, c, c, c};
marker = {'o', 's', '^', '*'};
for n = 1:Ndeg
    for m = 1:Nmesh
        pb = ParabolicBowl2d( order(n), ne(m), type );
        tic; pb.matSolve(); time(m, n) = toc;
        pos = makeNdgPostProcessFromNdgPhys( pb );
        tmp = pos.evaluateNormErr2(pb.fphys, pb.fext);
        err2(m, n) = tmp(1);
        tmp = pos.evaluateNormErr1(pb.fphys,pb.fext);
        err1(m, n) = tmp(1);
        tmp = pos.evaluateNormErrInf(pb.fphys, pb.fext);
        errInf(m, n) = tmp(1);

    end
    % print table
    fprintf('\n==================deg = %d==================\n', n);
    convergence_table(len, err1(:, n), err2(:, n), errInf(:, n), ...
        time(:, n))
    
    % plot figure
    co = color{n}; ma = marker{n};
    figure(1); plot(ne, err1(:, n), [co, ma, linestyle],...
        'LineWidth', linewidth, ...
        'MarkerSize', markersize ...
        ); 
    hold on;
    figure(2); plot(ne, err2(:, n), [co, ma, linestyle],...
        'LineWidth', linewidth, ...
        'MarkerSize', markersize ...
        ); 
    hold on;
    figure(3); plot(ne, errInf(:, n), [co, ma, linestyle],...
        'LineWidth', linewidth, ...
        'MarkerSize', markersize ...
        ); 
    hold on;
    
    figure(4); plot(dofs(:, n), time(:, n), [co, ma, linestyle],...
        'LineWidth', linewidth, ...
        'MarkerSize', markersize ...
        ); 
    hold on;
end

ylabel_str = {'$L_1$', '$L_2$', '$L_\infty$', '$time(s)$'};
fontsize = 16;
for n = 1:4
    figure(n);
    box on; grid on;
    set(gca, 'XScale', 'log', 'YScale', 'log');
    lendstr = cell(Ndeg, 1);
    for m = 1:Ndeg
        lendstr{m} = ['$p=', num2str(m), '$'];
    end
    legend(lendstr, 'box', 'off',...
        'Interpreter', 'Latex', 'FontSize', fontsize);
    xlabel('$DOFs$', 'Interpreter', 'Latex', 'FontSize', fontsize);
    ylabel(ylabel_str{n}, 'Interpreter', 'Latex', 'FontSize', fontsize);
end

end

function t1 = convergence_table(len, err1, err2, errInf, time)
t1 = table;
t1.len = len(:);
t1.('err1') = err1(:);
t1.('a1') = get_ratio(len, err1);

t1.('err2') = err2(:);
t1.('a2') = get_ratio(len, err2);

t1.('errf') = errInf(:);
t1.('af') = get_ratio(len, errInf);

t1.('time') = time(:);
end

function a = get_ratio(len, err)
Nmesh = numel(len);

a = zeros(Nmesh, 1);
for m = 2:Nmesh
    scal_ratio = log2( len(m)/len(m-1) );
    a(m) = log2( err(m)/err(m-1) )./scal_ratio;
end
end