function [rhsH, rhsQ] = SWERHS2d(mesh, h, q)

h = Q(:,:,1); hu = Q(:,:,2); hv = Q(:,:,3); g = 9.8;

%% flux
[Fx, Fy] = SWE2DFlux(Q);

%% numerical flux

hM  = zeros(Nfp(), mesh.nElement); hM(:) = h(mesh.vmapM);
hP  = zeros(Nfp(), mesh.nElement); hP(:) = h(mesh.vmapP); 
huM  = zeros(Nfp(), mesh.nElement); huM(:) = hu(mesh.vmapM); 
huP  = zeros(Nfp(), mesh.nElement); huP(:) = hu(mesh.vmapP);
hvM  = zeros(Nfp(), mesh.nElement); hvM(:) = hv(mesh.vmapM); 
hvP  = zeros(Nfp(), mesh.nElement); hvP(:) = hv(mesh.vmapP);

lamda1 = (mesh.nx.*huM./hM + mesh.ny.*hvM./hM) - sqrt(g.*hM);
lamda = lamda1;
lamda1 = (mesh.nx.*huP./hP + mesh.ny.*hvP./hP) - sqrt(g.*hP);
lamda(lamda1>lamda) = lamda1(lamda1>lamda);
% max{u}
lamda1 = (mesh.nx.*huM./hM + mesh.ny.*hvM./hM);
lamda(lamda1>lamda) = lamda1(lamda1>lamda);
lamda1 = (mesh.nx.*huP./hP + mesh.ny.*hvP./hP);
lamda(lamda1>lamda) = lamda1(lamda1>lamda);
% max{u+a}
lamda1 = (mesh.nx.*huM./hM + mesh.ny.*hvM./hM) + sqrt(g.*hM);
lamda(lamda1>lamda) = lamda1(lamda1>lamda);
lamda1 = (mesh.nx.*huP./hP + mesh.ny.*hvP./hP) + sqrt(g.*hP);
lamda(lamda1>lamda) = lamda1(lamda1>lamda);
%
lamda = abs(lamda);

QM(:,:,1) = hM; QM(:,:,2) = huM; QM(:,:,3)=hvM;
QP(:,:,1) = hP; QP(:,:,2) = huP; QP(:,:,3)=hvP;

[FxM, FyM] = SWE2DFlux(QM); [FxP,FyP] = SWE2DFlux(QP);

dF = zeros(Nfp(), mesh.nElement, 3);
dF(:,:,1) = mesh.nx.*(FxM(:,:,1) - FxP(:,:,1))./2 + mesh.ny.*(FyM(:,:,1) - FyP(:,:,1))./2 ...
    - lamda.*(QM(:,:,1) - QP(:,:,1))./2;
dF(:,:,2) = mesh.nx.*(FxM(:,:,2) - FxP(:,:,2))./2 + mesh.ny.*(FyM(:,:,2) - FyP(:,:,2))./2 ...
    - lamda.*(QM(:,:,2) - QP(:,:,2))./2;
dF(:,:,3) = mesh.nx.*(FxM(:,:,3) - FxP(:,:,3))./2 + mesh.ny.*(FyM(:,:,2) - FyP(:,:,2))./2 ...
    - lamda.*(QM(:,:,3) - QP(:,:,3))./2;

%% RHS
rhsQ(:,:,1) = -( Dr()*(mesh.rx.*Fx(:,:,1)+mesh.ry.*Fy(:,:,1))+ Ds()*(mesh.sx.*Fx(:,:,1) + mesh.sy.*Fy(:,:,1)) ) ...
    + invM()*(Fmat()*(mesh.fScale.*dF(:,:,1)));
rhsQ(:,:,2) = -( Dr()*(mesh.rx.*Fx(:,:,2)+mesh.ry.*Fy(:,:,2))+ Ds()*(mesh.sx.*Fx(:,:,2) + mesh.sy.*Fy(:,:,2)) ) ...
    + invM()*(Fmat()*(mesh.fScale.*dF(:,:,2)));
rhsQ(:,:,3) = -( Dr()*(mesh.rx.*Fx(:,:,3)+mesh.ry.*Fy(:,:,3))+ Ds()*(mesh.sx.*Fx(:,:,3) + mesh.sy.*Fy(:,:,3)) ) ...
    + invM()*(Fmat()*(mesh.fScale.*dF(:,:,3)));

end