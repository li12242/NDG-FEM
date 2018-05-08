%> @brief Initial function
%> This function will set the solver's mesh objects and initialize
%> the physical field
function initPhysFromOptions( obj, mesh )
% set mesh object
obj.meshUnion = mesh;
obj.Nmesh = numel(mesh);

% set option
obj.option = obj.setOption( obj.option );
% initilize physical field
obj.fphys = obj.setInitialField( );
end% func
