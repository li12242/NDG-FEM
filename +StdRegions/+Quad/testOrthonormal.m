% orithonormal test

% get values of orthonormal basis on interpolation points V(:, i)
np = 21;
nOrder = 4;
[r, s] = StdRegions.Quad.getCoor(np - 1);
V = StdRegions.Quad.getVandMatrix(4, r, s);

% get Gauss integral coefficients
[~,w] = Polylib.zwglj(np);
[w1, w2] = meshgrid(w, w);
w1 = w1(:); w2 = w2(:);

% get two polynomial from basis
v1 = 12; v2 = 12;
temp1 = V(:, v1); temp2 = V(:, v2);

% get the intergal result
sum(w1.*temp1.*w2.*temp2)