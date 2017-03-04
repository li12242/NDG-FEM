% orithonormal test

% get values of orthonormal basis on interpolation points V(:, i)
np = 5;
nOrder = 4;
[r, s] = StdRegions.Quad.getCoor(np - 1);
V = StdRegions.Quad.getVandMatrix(nOrder, r, s);

% get Gauss integral coefficients
[~,w] = Polylib.zwglj(np);
[w1, w2] = meshgrid(w, w);
w1 = w1(:); w2 = w2(:);

% get derivative matrix
[Dr,Ds,Drw,Dsw] = StdRegions.Quad.getDeriMatrix(nOrder,r,s,V);

Dr*r
Dr*s
Ds*r
Ds*s