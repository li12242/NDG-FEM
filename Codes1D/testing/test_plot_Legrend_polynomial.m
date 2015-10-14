N = 8; Np = N+1; 
% Nfp = 1; Nfaces=2;

% Compute basic Legendre Gauss Lobatto grid
% r = JacobiGL(0,0,N);

% Build reference element matrices
% V  = Vandermonde1D(N, r); invV = inv(V);

row = floor(sqrt(Np)); col = ceil(Np/row);

figure('Color', 'w')
for i = 1:Np
    x = linspace(-1, 1, 100)';
    P = JacobiP(x, 0, 0, i-1);
    subplot(row, col, i)
    plot(x, P)
    
    xlabel(['P_',num2str(i-1)])
end

