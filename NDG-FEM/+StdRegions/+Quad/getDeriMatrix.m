function [Dr,Ds,Drw,Dsw] = getDeriMatrix(N,r,s,V)
% function [Dr,Ds] = Dmatrices2D(N,r,s,V)
% Purpose : Initialize the (r,s) differentiation matrices on quadrilateral

%%
% $Dr_{ij} = \frac{\partial \varphi(\mathbf{\xi})}{\partial r}$
% $Ds_{ij} = \frac{\partial \varphi(\mathbf{\xi})}{\partial s}$

[Vr,Vs] = GradVandermonde2D(N,r,s);
Dr = Vr/V; Ds = Vs/V;

Drw = (V*Vr')/(V*V'); Dsw = (V*Vs')/(V*V');
end% func


function [V2Dr,V2Ds] = GradVandermonde2D(nOrder,r,s)

V2Dr = zeros(numel(r), (nOrder+1)^2 );
V2Ds = zeros(numel(r), (nOrder+1)^2 );
np = nOrder + 1; % number of points on edge
for i = 0:nOrder
    % P_{j-1}(r_i)$
    temp = Polylib.GradJacobiP(r(:), 0, 0, i);
    V2Dr(:,i*np+1:(i+1)*np) = repmat(temp, 1, np);
    for j = 0:nOrder
        V2Dr(:,i*np+1 + j) = V2Dr(:,i*np+1 + j).*Polylib.JacobiP(s(:), 0, 0, j);
    end% for
end% for

for i = 0:nOrder
    % P_{j-1}(r_i)$
    temp = Polylib.JacobiP(r(:), 0, 0, i);
    V2Ds(:,i*np+1:(i+1)*np) = repmat(temp, 1, np);
    for j = 0:nOrder
        V2Ds(:,i*np+1 + j) = V2Ds(:,i*np+1 + j).*Polylib.GradJacobiP(s(:), 0, 0, j);
    end% for
end% for
end% func