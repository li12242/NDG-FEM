function [c, q] = StyDiffInit(mesh)
% 1D Steady diffusion
% Purpose: initial condition
c = zeros(size(mesh.x));
q = zeros(size(mesh.x));
end% func