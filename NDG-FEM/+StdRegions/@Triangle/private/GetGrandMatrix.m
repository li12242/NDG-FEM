function [ Dr, Ds ] = GetGrandMatrix( N,r,s,V )
%GETGRANDMATRIX ����ڵ�������������꣨r,s���ĵ�������
%   ���е����������� 
%   Dr_{ij} = \frac{\partial l_j}{\partial r} |_{(r = r_i, s = s_i)}

Np = (N+1)*(N+2)/2;
Vdr = zeros(Np, Np);
Vds = zeros(Np, Np);
for ind = 1:Np
    [Vdr(:,ind), Vds(:,ind)] = GradOrthogonalFun( r,s,N,ind );
end% for

Dr = Vdr/V; Ds = Vds/V;
end
