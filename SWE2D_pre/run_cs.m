function run_cs(type)
%RUN_CS Summary of this function goes here
%   Detailed explanation goes here

% check cell type
if (type ~= ndg_lib.std_cell_type.Tri) ...
        && (type ~= ndg_lib.std_cell_type.Quad)
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
cs_ext = load('SWE2d/@swe2d_cs/cs_ext.mat');

errInf = zeros(Nmesh, Ndeg, Nfield);
err2 = zeros(Nmesh, Ndeg, Nfield);
err1 = zeros(Nmesh, Ndeg, Nfield);
quad_type = ndg_lib.std_cell_type.Quad;
linewidth = 1.5; 
markersize = 8;
color = {'b', 'r', 'g', 'm'};
marker = {'o', 's', '^', '*'};
linestyle = '--';

for n = 1:Ndeg
    for m = 1:Nmesh
        cs = swe2d_cs(order(n), ne(m), quad_type);
        cs.init; 
        tic; cs.RK45; time(m, n) = toc;
        [err1(m, n, :), err2(m, n, :), errInf(m, n, :)]...
            = norm_err(cs, cs_ext.cs);
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
    
    figure(4); plot(ne, time(:, n), [co, ma, linestyle],...
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

function [err1, err2, errInf] = norm_err(cs, cs_ext)
% detector = ndg_utility.detector.detector2d(cs.mesh, ...
%     cs_ext.mesh.x(:)', cs_ext.mesh.y(:)', cs.ftime, cs.ftime, cs.Nfield);
% detector.collect( cs.f_Q, cs.ftime );
% f_ext(:,:,1) = reshape(detector.dQ(:,:,1), cs_ext.mesh.cell.Np, cs_ext.mesh.K);
% f_ext(:,:,2) = reshape(detector.dQ(:,:,2), cs_ext.mesh.cell.Np, cs_ext.mesh.K);
% f_ext(:,:,3) = reshape(detector.dQ(:,:,3), cs_ext.mesh.cell.Np, cs_ext.mesh.K);
detector = ndg_utility.detector.detector2d(cs_ext.mesh, ...
    cs.mesh.x(:)', cs.mesh.y(:)', cs.ftime, cs.ftime, cs.Nfield);
detector.collect( cs_ext.f_Q, cs.ftime );
f_ext(:,:,1) = reshape(detector.dQ(:,:,1), cs.mesh.cell.Np, cs.mesh.K);
f_ext(:,:,2) = reshape(detector.dQ(:,:,2), cs.mesh.cell.Np, cs.mesh.K);
f_ext(:,:,3) = reshape(detector.dQ(:,:,3), cs.mesh.cell.Np, cs.mesh.K);


f_abs = cs.f_Q - f_ext;
area = sum(cs.mesh.vol);

err1 = zeros(cs.Nfield, 1);
err2 = zeros(cs.Nfield, 1);
errInf = zeros(cs.Nfield, 1);

for fld = 1:cs.Nfield
    temp = abs( f_abs(:,:,fld) );
    err1(fld) = sum( ...
        cs.mesh.cell_mean(temp).*cs.mesh.vol )./area;
    errInf(fld) = max( max(temp) );
    
    temp = f_abs(:,:,fld).*f_abs(:,:,fld);
    err2(fld) = sqrt( sum( ...
        cs.mesh.cell_mean(temp).*cs.mesh.vol ) )./area;
end
end% func

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

