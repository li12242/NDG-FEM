classdef NdgAbstractAdvSubSolver < handle
    
    properties( SetAccess = protected )
        regionId
        phys
        Nmesh
    end
    
    properties
        % length/area/volume of each control volume
        LAVToCV
        % outword normal vector of inner edge
        nxInner
        % outword normal vector of inner edge
        nyInner
        % outword normal vector of inner edge
        nzInner
        % length/area of inner edge
        flen
        % local node index of each inner edge
        n1
        % adjacent node index of each inner edge
        n2
        % outword normal vector of facial edge
        nx
        % outword normal vector of facial edge
        ny
        % outword normal vector of facial edge
        nz
        % length/area of facial edge
        Js
    end
    
    methods
        function obj = NdgAbstractAdvSubSolver( phys, regionId )
            obj.phys = phys;
            obj.regionId = regionId;
            obj.Nmesh = phys.Nmesh;
            
            obj.LAVToCV = cell( obj.Nmesh, 1 );
            obj.nxInner = cell( obj.Nmesh, 1 );
            obj.nyInner = cell( obj.Nmesh, 1 );
            obj.nzInner = cell( obj.Nmesh, 1 );
            obj.flen = cell( obj.Nmesh, 1 );
            obj.nx = cell( obj.Nmesh, 1 );
            obj.ny = cell( obj.Nmesh, 1 );
            obj.nz = cell( obj.Nmesh, 1 );
            obj.Js = cell( obj.Nmesh, 1 );
            
            obj.n1 = cell( obj.Nmesh, 1 );
            obj.n2 = cell( obj.Nmesh, 1 );
            for m = 1:obj.Nmesh
                mesh = phys.meshUnion(m);
                obj.n1{m} = mesh.cell.n1;
                obj.n2{m} = mesh.cell.n2;
                
                LAVToCV = mesh.cell.LAVToCV...
                    ./phys.meshUnion(m).cell.LAV;
                obj.LAVToCV{m} = bsxfun(@times, mesh.LAV, LAVToCV);
                
                [ obj.nxInner{m}, obj.nyInner{m}, obj.nzInner{m}, obj.flen{m} ] = ...
                    mesh.cell.matEvaluateInnerControlEdgeInfo( mesh.x, mesh.y, mesh.z );
                
                [ obj.nx{m}, obj.ny{m}, obj.nz{m}, obj.Js{m} ] = ...
                    mesh.cell.assembleNormalVector( mesh.x, mesh.y, mesh.z );
                
            end
        end
    end
    
    methods( Abstract )
        evaluateAdvectionRHS( obj, fphys )
    end
    
end

