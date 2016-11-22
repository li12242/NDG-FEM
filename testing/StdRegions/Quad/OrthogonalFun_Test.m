tol = 1e-10;

Nmax = 10; Np = Nmax+1;
[z,w] = Polylib.zwglj(Nmax+1);

% 权重
temp = w*w';
wei = temp(:);
J = 4/sum(wei);

% 高斯积分节点
r = z*ones(1, Np);
s = ones(Np, 1)*z';
r = r(:);
s = s(:);

%% Quad: N = 3
N= 3;
quad = StdRegions.Quad(N);
V = quad.GetVandMatrix(N, r, s);

for i = 1:quad.nNode
    for j = 1:quad.nNode
        temp = J*sum(wei.*V(:,i).*V(:,j));
        assert( abs(temp - (i==j))<tol )
    end
end

%% Quad: N = 4
N = 4;
quad = StdRegions.Quad(N);
V = quad.GetVandMatrix(N, r, s);

for i = 1:quad.nNode
    for j = 1:quad.nNode
        temp = J*sum(wei.*V(:,i).*V(:,j));
        assert( abs(temp - (i==j))<tol )
    end
end

%% Quad: N = 5
N = 5;
quad = StdRegions.Quad(N);
V = quad.GetVandMatrix(N, r, s);

for i = 1:quad.nNode
    for j = 1:quad.nNode
        temp = J*sum(wei.*V(:,i).*V(:,j));
        assert( abs(temp - (i==j))<tol )
    end
end