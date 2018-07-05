function Fmask = AssembleFacialNodeIndex(obj)

maxnfp = max(obj.Nfp);
Fmask = zeros(maxnfp, obj.Nface);

% vertical faces
f = 1;
ind = find( abs(obj.s + 1) < 1e-10 );
Fmask(1:obj.Nfp(f), f) = ind;

f = 2;
ind = find( abs(obj.r + obj.s) < 1e-10 );
Fmask(1:obj.Nfp(f), f) = ind;

f = 3;
ind = find( abs(obj.r + 1) < 1e-10 );
Fmask(1:obj.Nfp(f), f) = ind;
% horizontal faces
f = 4;
ind = find( abs(obj.t + 1) < 1e-10 );
Fmask(1:obj.Nfp(f), f) = ind;

f = 5;
ind = find( abs(obj.t - 1) < 1e-10 );
Fmask(1:obj.Nfp(f), f) = ind;

end% funca