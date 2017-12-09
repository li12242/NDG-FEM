classdef NdgQuadFreeStrongFormSolver < handle
    %NDGQUADFREESTRONGFORMSOLVER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties( SetAccess = protected )
        Dr
        Ds
        Dt
        LIFT
        rx
        ry
        rz
        sx
        sy
        sz
        tx
        ty
        tz
        J
        nx
        ny
        nz
        Js
    end
    
    methods
        function obj = NdgQuadFreeStrongFormSolver( phys )
            
            obj.Dr = cell( phys.Nmesh, 1 );
            obj.Ds = cell( phys.Nmesh, 1 );
            obj.Dt = cell( phys.Nmesh, 1 );
            obj.LIFT = cell( phys.Nmesh, 1 );
            obj.rx = cell( phys.Nmesh, 1 );
            obj.ry = cell( phys.Nmesh, 1 );
            obj.rz = cell( phys.Nmesh, 1 );
            obj.sx = cell( phys.Nmesh, 1 );
            obj.sy = cell( phys.Nmesh, 1 );
            obj.sz = cell( phys.Nmesh, 1 );
            obj.tx = cell( phys.Nmesh, 1 );
            obj.ty = cell( phys.Nmesh, 1 );
            obj.tz = cell( phys.Nmesh, 1 );
            obj.J = cell( phys.Nmesh, 1 );
            obj.nx = cell( phys.Nmesh, 1 );
            obj.ny = cell( phys.Nmesh, 1 );
            obj.nz = cell( phys.Nmesh, 1 );
            obj.Js = cell( phys.Nmesh, 1 );
            
            for m = 1:phys.Nmesh
                mesh = phys.meshUnion( m );
                [ obj.Dr{m}, obj.Ds{m}, obj.Dt{m} ] = obj.assembleDerivativeMatrix( mesh.cell );
                [ obj.LIFT{m} ] = obj.assembleLiftMatrix( mesh.cell );
                [ obj.rx{m}, obj.ry{m}, obj.rz{m}, ...
                    obj.sx{m}, obj.sy{m}, obj.sz{m}, ...
                    obj.tx{m}, obj.ty{m}, obj.tz{m}, obj.J{m} ] ...
                    = obj.assembleJacobianFactor( mesh );
                [ obj.nx{m}, obj.ny{m}, obj.nz{m}, obj.Js{m} ] = obj.assembleNormalVector( mesh );
            end
        end
    end
    
    methods( Static, Access = private )
        function [ rx, ry, rz, sx, sy, sz, tx, ty, tz, J ] = assembleJacobianFactor( mesh )
            [ rx, ry, rz, sx, sy, sz, tx, ty, tz, J ] ...
                = mesh.cell.assembleJacobianMatrix( mesh.x, mesh.y, mesh.z );
        end
        
        function [ nx, ny, nz, Js ] = assembleNormalVector( mesh )
            [ nx, ny, nz, Js ] = mesh.cell.assembleNormalVector( mesh.x, mesh.y, mesh.z );
        end
        
        function [ Dr, Ds, Dt ] = assembleDerivativeMatrix( cell )
            Dr = cell.Dr;
            Ds = cell.Ds;
            Dt = cell.Dt;
        end
        
        function LIFT = assembleLiftMatrix( cell )
            Mes = zeros(cell.Np, cell.TNfp);
            sk = 1;
            for f = 1:cell.Nface
                fcell = getStdCell(cell.N, cell.faceType(f));
                row = cell.Fmask(:, f);
                row = row(row ~= 0);
                for n = 1:fcell.Np
                    Mes(row, sk) = fcell.M(:, n);
                    sk = sk + 1;
                end
            end
            LIFT = cell.invM * Mes;
        end
    end
    
end

