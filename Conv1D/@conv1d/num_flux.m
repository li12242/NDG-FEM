function flux = num_flux(obj, cM, uM, cP, uP, nx)
%NUM_FLUX Summary of this function goes here
%   Detailed explanation goes here

flux = zeros(obj.mesh.cell.Nfptotal, obj.mesh.K);
ind = (nx.*uM >= 0); flux(ind) = nx(ind).*uM(ind).*cM(ind);
ind = (nx.*uM <  0); flux(ind) = nx(ind).*uP(ind).*cP(ind);

end

