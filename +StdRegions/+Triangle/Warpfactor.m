function warp = Warpfactor(N, rout)
% Compute scaled warp function at order N based on rout interpolation nodes

% Compute LGL and equidistant node distribution
[LGLr,~] = Polylib.zwglj(N+1); req = linspace(-1,1,N+1)';

% Compute V based on req
%Veq = StdRegions.Line.GetVandMatrix(N, req);
obj=StdRegions.Line(N);
Veq = obj.GetVandMatrix(N, req);
% Evaluate Lagrange polynomial at rout
Nr = length(rout); %Lmat=zeros(N+1,Nr); 
Pmat = zeros(N+1,Nr);
for i=1:N+1
  Pmat(i,:) = Polylib.JacobiP(rout, 0, 0, i-1);
end
Lmat = Veq'\Pmat;
% Compute warp factor
warp = Lmat'*(LGLr - req);
% Scale factor
zerof = (abs(rout)<1.0-1.0e-10); sf = 1.0 - (zerof.*rout).^2;
warp = warp./sf;
end