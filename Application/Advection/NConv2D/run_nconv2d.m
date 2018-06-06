function run_nconv2d(type)
%RUN_NCONV2D Summary of this function goes here
%   Detailed explanation goes here

% check cell type
if (type == ndg_lib.std_cell_type.Quad)
    casename = 'NConv2D/@nconv2d_sin/mesh/quad500/mesh';
elseif (type == ndg_lib.std_cell_type.Tri)
    casename = 'NConv2D/@nconv2d_sin/mesh/tri1000/mesh';
else
    error(['Illegal cell type: ', num2str(type)]);
end
% parameter
order = [1,2,3];
ne = [20, 40, 80];
len = 20./ne;

Nmesh = numel(ne);
Ndeg = numel(order);
Nfield = 3;
dofs = zeros(Nmesh, Ndeg);
time = zeros(Nmesh, Ndeg);
for n = 1:Ndeg
    for m = 1:Nmesh
        if (type == ndg_lib.std_cell_type.Tri)
            dofs(m, n) = ne(m).^2 *2 * (order(n)+1)*(order(n)+2)/2 ;
        elseif (type ~= ndg_lib.std_cell_type.Tri)
            dofs(m, n) = ne(m).^2 * (order(n)+1).^2;
        end
    end
end
% get the result on fine mesh as the exact solution

errInf = zeros(Nmesh, Ndeg, Nfield);
err2 = zeros(Nmesh, Ndeg, Nfield);
err1 = zeros(Nmesh, Ndeg, Nfield);
linewidth = 1.5; 
markersize = 8;
% color = {'b', 'r', 'g', 'm'};
color = {'k', 'k', 'k', 'k'};
marker = {'o', 's', '^', '*'};
linestyle = '-';

for n = 1:Ndeg
    for m = 1:Nmesh
        if m == 1
            nconv = nconv2d_sin(order(n), casename, type);
        else
            nconv.refine_mesh;
        end
%         nconv = nconv2d_sin(order(n), ne(m), type);
        nconv.init; 
        tic; nconv.RK45; time(m, n) = toc;
        err2(m, n) = nconv.norm_err2(nconv.ftime);
        err1(m, n) = nconv.norm_err1(nconv.ftime);
        errInf(m, n) = nconv.norm_errInf(nconv.ftime);
    end
    % print table
    fprintf('\n==================deg = %d==================\n', n);
    t = convergence_table(len, err1(:, n, 1), err2(:, n, 1), errInf(:, n, 1), ...
        time(:, n))
    
    % plot figure
    co = color{n}; ma = marker{n};
    figure(1); plot(dofs(:, n), err1(:, n, 1), [co, ma, linestyle],...
        'LineWidth', linewidth, ...
        'MarkerSize', markersize ...
        ); 
    hold on;
    figure(2); plot(dofs(:, n), err2(:, n, 1), [co, ma, linestyle],...
        'LineWidth', linewidth, ...
        'MarkerSize', markersize ...
        ); 
    hold on;
    figure(3); plot(dofs(:, n), errInf(:, n, 1), [co, ma, linestyle],...
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
