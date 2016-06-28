function [ux,uy] = Grad2D(mesh, u)
%% Compute 2D gradient field of scalar field U
% Where $\nabla u=\left( \frac{\partial u}{\partial x}, \frac{\partial u}{\partial y} \right)$ 

shape = mesh.Shape;

ur = shape.Dr*u; 
us = shape.Ds*u;

ux = mesh.rx.*ur + mesh.sx.*us; 
uy = mesh.ry.*ur + mesh.sy.*us;
return
