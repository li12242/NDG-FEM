function [ Dr,Ds ] = GetGrandMatrix( N,r,s,V )
%GETGRANDMATRIX Summary of this function goes here
%   Detailed explanation goes here

Np = (N+1)*(N+1);
Vdr = zeros(Np, Np);
Vds = zeros(Np, Np);
for ind = 1:Np
    [Vdr(:,ind), Vds(:,ind)] = GradOrthogonalFun( r,s,N,ind );
end% for

Dr = Vdr/V; Ds = Vds/V;
end

