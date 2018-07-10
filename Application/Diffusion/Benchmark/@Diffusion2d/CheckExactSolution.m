function CheckExactSolution( obj )
%CHE Summary of this function goes here
%   Detailed explanation goes here

mesh = obj.meshUnion;
visual = Visual2d( mesh );
fext = obj.getExtFunc( mesh.x, mesh.y, obj.getOption('finalTime') );
delta = obj.fphys{1} - fext;

visual.drawResult( delta );
end

