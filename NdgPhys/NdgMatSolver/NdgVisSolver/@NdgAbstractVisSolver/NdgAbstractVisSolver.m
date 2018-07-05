classdef NdgAbstractVisSolver < handle
    
    properties( SetAccess = protected )
        % num of mesh
        Nmesh
        % physical object
        phys
        % auxiliary variabel
        px, py, pz
        % variable index, and its corresponding rhs index
        varId, rhsId
        % num of field
        Nfield
    end
    
    properties( SetAccess = public )
        % viscosity
        mx, my, mz
    end
    
    methods
        function obj = NdgAbstractVisSolver( phys, varId, rhsId )
            obj.phys = phys;
            obj.varId = varId;
            obj.rhsId = rhsId;
            obj.Nfield = numel( varId );
            
            obj.Nmesh = phys.Nmesh;
            obj.px = cell( obj.Nmesh, 1 );
            obj.py = cell( obj.Nmesh, 1 );
            obj.pz = cell( obj.Nmesh, 1 );
            
            obj.mx = cell( obj.Nmesh, 1 );
            obj.my = cell( obj.Nmesh, 1 );
            obj.mz = cell( obj.Nmesh, 1 );
            for m = 1:obj.Nmesh
                obj.px{m} = zeros( phys.meshUnion(m).cell.Np, phys.meshUnion(m).K );
                obj.py{m} = zeros( phys.meshUnion(m).cell.Np, phys.meshUnion(m).K );
                obj.pz{m} = zeros( phys.meshUnion(m).cell.Np, phys.meshUnion(m).K );
            end
        end% func
    end% methods
    
    methods(Abstract)
        matEvaluateRHS( obj, fphys, frhs );
    end% methods
    
end
