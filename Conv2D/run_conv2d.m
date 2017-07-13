function run_conv2d
casename{1} = 'Conv2D/@conv2d_diffusion/mesh/quad_500/quad_500';
K = 505;
n = [0, 1, 2];
ne = K*4.^n;
len = 0.5.^n;
k = [1, 2, 3];
% ne = [20, 40, 60, 80];
% len = 1./ne;

Nmesh = numel(ne);
Ndeg = numel(k);
dofs = zeros(Nmesh, Ndeg);
time = zeros(Nmesh, Ndeg);
for n = 1:Ndeg
    for m = 1:Nmesh
        dofs(m, n) = ne(m) * (k(n)+1).^2;
    end
end

errInf = zeros(Nmesh, Ndeg);
err2 = zeros(Nmesh, Ndeg);
err1 = zeros(Nmesh, Ndeg);
quad_type = ndg_lib.std_cell_type.Quad;
linewidth = 1.5; 
markersize = 8;
color = {'b', 'r', 'g', 'm'};
marker = {'o', 's', '^', '*'};
linestyle = '--';
for n = 1:Ndeg
    for m = 1:Nmesh
        if (m == 1)
            conv = conv2d_diffusion(k(n), casename{m}, quad_type);
        else
            conv.refine_mesh(1);
        end
        %conv.mesh.J = repmat(mean(conv.mesh.J), conv.mesh.cell.Np, 1);
        conv.init; 
        tic; conv.RK45_solve; time(m, n) = toc;
        err2(m, n) = conv.norm_err2(conv.ftime);
        err1(m, n) = conv.norm_err1(conv.ftime);
        errInf(m, n) = conv.norm_errInf(conv.ftime);
    end
    % print table
    fprintf('\n==================deg = %d==================\n', n);
    t = convergence_table(len, err1(:, n), err2(:, n), errInf(:, n), ...
        time(:, n))
    
    % plot figure
    co = color{n}; ma = marker{n};
    figure(1); plot(dofs(:, n), err1(:, n), [co, ma, linestyle],...
        'LineWidth', linewidth, ...
        'MarkerSize', markersize, ...
        'MarkerFaceColor', co); 
    hold on;
    figure(2); plot(dofs(:, n), err2(:, n), [co, ma, linestyle],...
        'LineWidth', linewidth, ...
        'MarkerSize', markersize, ...
        'MarkerFaceColor', co); 
    hold on;
    figure(3); plot(dofs(:, n), errInf(:, n), [co, ma, linestyle],...
        'LineWidth', linewidth, ...
        'MarkerSize', markersize, ...
        'MarkerFaceColor', co); 
    hold on;
    
    figure(4); plot(dofs(:, n), time(:, n), [co, ma, linestyle],...
        'LineWidth', linewidth, ...
        'MarkerSize', markersize, ...
        'MarkerFaceColor', co); 
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