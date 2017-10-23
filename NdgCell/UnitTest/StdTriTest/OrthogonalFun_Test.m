% 验证正交基函数
tol = 1e-10;

% 三角形数值积分公式参考自
% <<http://people.sc.fsu.edu/~jburkardt/datasets/
% quadrature_rules_tri/quadrature_rules_tri.html>>

% 顶点坐标
vert = load('GaussQaud/triGauss8x8_r.txt');
vx = vert(:, 1);
vy = vert(:, 2);

% 高斯积分点坐标
points = load('GaussQaud/triGauss8x8_x.txt');
x = points(:, 1); 
y = points(:, 2);

% 权重系数
w = load('GaussQaud/triGauss8x8_w.txt');
J = 2/sum(w); % Jacobi 转换系数 - 标准三角形单元面积为2

% 转换为面积坐标
L1 = y;
L3 = x;
L2 = 1-L1-L3;

% 根据面积坐标计算对应标准单元内节点坐标
r = -L2 + L3 - L1; s = -L2 - L3 + L1;

% 最大积分阶数
MaxDeg = 64;

%% Tri: N = 2
N = 2;
tri = StdRegions.Triangle(N);
V = tri.GetVandMatrix(N, r, s);

for i = 1:tri.nNode
    for j = 1:tri.nNode
        temp = J*sum(w.*V(:,i).*V(:,j));
        assert( abs(temp - (i==j))<tol )
    end
end

%% Tri: N = 3
N = 3;
tri = StdRegions.Triangle(N);
V = tri.GetVandMatrix(N, r, s);

for i = 1:tri.nNode
    for j = 1:tri.nNode
        temp = J*sum(w.*V(:,i).*V(:,j));
        assert( abs(temp - (i==j))<tol )
    end
end

%% Tri: N = 4
N = 4;
tri = StdRegions.Triangle(N);
V = tri.GetVandMatrix(N, r, s);

for i = 1:tri.nNode
    for j = 1:tri.nNode
        temp = J*sum(w.*V(:,i).*V(:,j));
        assert( abs(temp - (i==j))<tol )
    end
end

%% Tri: N = 5
N = 5;
tri = StdRegions.Triangle(N);
V = tri.GetVandMatrix(N, r, s);

for i = 1:tri.nNode
    for j = 1:tri.nNode
        temp = J*sum(w.*V(:,i).*V(:,j));
        assert( abs(temp - (i==j))<tol )
    end
end

%% Tri: N = 6
N = 6;
tri = StdRegions.Triangle(N);
V = tri.GetVandMatrix(N, r, s);

for i = 1:tri.nNode
    for j = 1:tri.nNode
        temp = J*sum(w.*V(:,i).*V(:,j));
        assert( abs(temp - (i==j))<tol )
    end
end

%% Tri: N = 7
N = 7;
tri = StdRegions.Triangle(N);
V = tri.GetVandMatrix(N, r, s);

for i = 1:tri.nNode
    for j = 1:tri.nNode
        temp = J*sum(w.*V(:,i).*V(:,j));
        assert( abs(temp - (i==j))<tol )
    end
end