function [rhsQ, lamda] = SWE1DRHS_LF(mesh, Q, Nfp, Dr, Fmat, invM, dt)
% 

% Interface function handles, 目的减少对象调用时超长API
% Nfp = @(x)mesh.Shape.nFaceNode;
% Dr = @(x)mesh.Shape.Dr;
% Fmat = @(x)mesh.Shape.FaceMassMatrixSmall;
% invM = @(x)mesh.Shape.invM;

h = Q(:,:,1); hu = Q(:,:,2);
g = 9.8;

%% flux
F = SWEFlux(Q);

%% numerical flux

hM  = zeros(Nfp(), mesh.nElement); hM(:) = h(mesh.vmapM);
hP  = zeros(Nfp(), mesh.nElement); hP(:) = h(mesh.vmapP); 
huM  = zeros(Nfp(), mesh.nElement); huM(:) = hu(mesh.vmapM); 
huP  = zeros(Nfp(), mesh.nElement); huP(:) = hu(mesh.vmapP);

% imply boundary condiation
% dudx = mesh.rx(:,end).*( Dr()*(hu(:,end)./h(:,end)) );
% dudx = dudx(end);
% u1 = (hu(end,end)./h(end,end)) ...
%     - dt*( (hu(end,end)./h(end,end) + sqrt(g*h(end,end))) )*dudx;
% dhdx = mesh.rx(:,end).*( Dr()*h(:,end) );
% dhdx = dhdx(end);
% h1 = h(end,end) - dt*( (hu(end,end)./h(end,end) + sqrt(g*h(end,end))) )*dhdx;
% hP(mesh.mapO) = h1;
% huP(mesh.mapO) = hP(mesh.mapO)*u1;
hP(mesh.mapO) = 2;
huP(mesh.mapO) = hP(mesh.mapO).*(huM(mesh.mapO)./hM(mesh.mapO) + ...
    2*sqrt(g*hM(mesh.mapO)) - 2*sqrt(g*hP(mesh.mapO)) );
% huP(mesh.mapO) = 0;

% lamda = max{u-a, u, u+a}, a=sqrt(g*h)
% lamda = zeros(Nfp(), mesh.nElement);
% max{u-a}
lamda1 = (huM./hM) - sqrt(g.*hM);
lamda = lamda1;
lamda1 = (huP./hP) - sqrt(g.*hP);
lamda(lamda1>lamda) = lamda1(lamda1>lamda);
% max{u}
lamda1 = (huM./hM);
lamda(lamda1>lamda) = lamda1(lamda1>lamda);
lamda1 = (huP./hP);
lamda(lamda1>lamda) = lamda1(lamda1>lamda);
% max{u+a}
lamda1 = (huM./hM) + sqrt(g.*hM);
lamda(lamda1>lamda) = lamda1(lamda1>lamda);
lamda1 = (huP./hP) + sqrt(g.*hP);
lamda(lamda1>lamda) = lamda1(lamda1>lamda);
%
lamda = abs(lamda);

QM(:,:,1) = hM; QM(:,:,2) = huM;
QP(:,:,1) = hP; QP(:,:,2) = huP;

FM = SWEFlux(QM); FP = SWEFlux(QP);

dF = zeros(Nfp(), mesh.nElement, 2);

dF(:,:,1) = mesh.nx.*(FM(:,:,1) - FP(:,:,1))./2 - lamda.*(QM(:,:,1) - QP(:,:,1))./2;
dF(:,:,2) = mesh.nx.*(FM(:,:,2) - FP(:,:,2))./2 - lamda.*(QM(:,:,2) - QP(:,:,2))./2;

%% RHS

rhsQ(:,:,1) = -mesh.rx.*(Dr()*F(:,:,1)) + invM()*(Fmat()*(mesh.fScale.*dF(:,:,1)));
rhsQ(:,:,2) = -mesh.rx.*(Dr()*F(:,:,2)) + invM()*(Fmat()*(mesh.fScale.*dF(:,:,2)));

end