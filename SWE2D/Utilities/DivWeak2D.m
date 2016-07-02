function [divu] = DivWeak2D(mesh, u,v)
%% Compute the 2D divergence of the vectorfield (u,v)
% Where $\nabla\cdot\mathbf{U}=\frac{\partial u}{\partial x}+\frac{\partial v}{\partial y}$ 

shape = mesh.Shape;

ur = shape.Drw*u; us = shape.Dsw*u; 
vr = shape.Drw*v; vs = shape.Dsw*v;

divu = mesh.rx.*ur + mesh.sx.*us + mesh.ry.*vr + mesh.sy.*vs;
return;
