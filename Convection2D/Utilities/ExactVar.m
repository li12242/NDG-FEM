function ev = ExactVar(x, y, t)

sigma = 125*1e3/33^2; 
xc = 0; yc = 3/5;
ev = exp(-sigma.*( (x - xc).^2 + (y - yc).^2) );

end% func