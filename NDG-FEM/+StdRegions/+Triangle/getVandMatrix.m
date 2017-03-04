function VandMatrix = getVandMatrix(N, r, s)
VandMatrix = zeros(length(r),(N+1)*(N+2)/2);
% Transfer to (a,b) coordinates
[a, b] = StdRegions.Triangle.rstoab(r, s);
% build the Vandermonde matrix
sk = 1;
for i=0:N
  for j=0:N - i
    VandMatrix(:,sk) = StdRegions.Triangle.Simplex2DP(a,b,i,j);
    sk = sk+1;
  end
end
end