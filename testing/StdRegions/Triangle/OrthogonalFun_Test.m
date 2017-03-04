% ��֤����������
tol = 1e-10;

% ��������ֵ���ֹ�ʽ�ο���
% <<http://people.sc.fsu.edu/~jburkardt/datasets/
% quadrature_rules_tri/quadrature_rules_tri.html>>

% ��������
vert = load('TriGaussQaud/triGauss8x8_r.txt');
vx = vert(:, 1);
vy = vert(:, 2);

% ��˹���ֵ�����
points = load('TriGaussQaud/triGauss8x8_x.txt');
x = points(:, 1); 
y = points(:, 2);

% Ȩ��ϵ��
w = load('TriGaussQaud/triGauss8x8_w.txt');
J = 2/sum(w); % Jacobi ת��ϵ�� - ��׼�����ε�Ԫ���Ϊ2

% ת��Ϊ�������
L1 = y;
L3 = x;
L2 = 1-L1-L3;

% ���������������Ӧ��׼��Ԫ�ڽڵ�����
r = -L2 + L3 - L1; s = -L2 - L3 + L1;

% �����ֽ���
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