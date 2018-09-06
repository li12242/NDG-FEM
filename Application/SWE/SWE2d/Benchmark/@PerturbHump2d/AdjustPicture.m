path = '/Users/mac/BaiduYun/Article/A novel well-balanced quadrature-free nodal discontinuous Galerkin scheme for the shallow water equations/manuscript/figure/Perturbation/big perturbation/';
filename = {'PerHump_1', ...
    'PerHump_2', 'PerHump_3', ...
    'PerHump_4', 'PerHump_5'};

for n = 1:5
    open( [path, filename{n}, '.fig'] );
    axis equal; 
    set( gcf, 'Position', [219   553   560   252]);
    set( gca, 'CLim', [0.996, 1.008]);
    print( [path, filename{n}], '-depsc', '-r300' );
    xlabel('$x$ (m)', 'FontSize', 20, 'Interpreter', 'Latex');
    ylabel('$y$ (m)', 'FontSize', 20, 'Interpreter', 'Latex');
end

path = '/Users/mac/BaiduYun/Article/A novel well-balanced quadrature-free nodal discontinuous Galerkin scheme for the shallow water equations/manuscript/figure/Perturbation/small perturbation/';
filename = {'PerHump_1', ...
    'PerHump_2', 'PerHump_3', ...
    'PerHump_4', 'PerHump_5'};

for n = 1:5
    open( [path, filename{n}, '.fig'] );
    axis equal; 
    set( gcf, 'Position', [219   553   560   252]);
    set( gca, 'CLim', [0.9996, 1.0008]);
    print( [path, filename{n}], '-depsc', '-r300' );
    xlabel('$x$ (m)', 'FontSize', 20, 'Interpreter', 'Latex');
    ylabel('$y$ (m)', 'FontSize', 20, 'Interpreter', 'Latex');
end