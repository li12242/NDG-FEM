classdef NdgQuadFreeStrongFormSolver < handle
    
    properties ( SetAccess = protected )
        %> differential matrix
        Dr, Ds, Dt
        %> Jacobian matrix entries
        rx, ry, rz
        sx, sy, sz
        tx, ty, tz
        %> determination of Jacobian matrix
        J
    end
    
    methods
        function obj = NdgQuadFreeStrongFormSolver( phys )
            cell_struct = cell( phys.Nmesh, 1 );
            obj.Dr = cell_struct;
            obj.Ds = cell_struct;
            obj.Dt = cell_struct;
            % obj.LIFT = cell_struct;
            obj.rx = cell_struct;
            obj.ry = cell_struct;
            obj.rz = cell_struct;
            obj.sx = cell_struct;
            obj.sy = cell_struct;
            obj.sz = cell_struct;
            obj.tx = cell_struct;
            obj.ty = cell_struct;
            obj.tz = cell_struct;
            obj.J = cell_struct;
            % obj.nx = cell_struct;
            % obj.ny = cell_struct;
            % obj.nz = cell_struct;
            % obj.Js = cell_struct;
            
            for m = 1:phys.Nmesh
                mesh = phys.meshUnion( m );
                obj.Dr{m} = mesh.cell.Dr;
                obj.Ds{m} = mesh.cell.Ds;
                obj.Dt{m} = mesh.cell.Dt;

                obj.rx{m} = mesh.rx; obj.ry{m} = mesh.ry; obj.rz{m} = mesh.rz;
                obj.sx{m} = mesh.sx; obj.sy{m} = mesh.sy; obj.sz{m} = mesh.sz;
                obj.tx{m} = mesh.tx; obj.ty{m} = mesh.ty; obj.tz{m} = mesh.tz;
                obj.J{m} = mesh.J;
                % [ obj.rx{m}, obj.ry{m}, obj.rz{m}, ...
                %     obj.sx{m}, obj.sy{m}, obj.sz{m}, ...
                %     obj.tx{m}, obj.ty{m}, obj.tz{m}, obj.J{m} ] ...
                %     = obj.assembleJacobianFactor( mesh );
            end
        end
    end
    
    % methods( Static, Access = private )
    %     function [ rx, ry, rz, sx, sy, sz, tx, ty, tz, J ] = assembleJacobianFactor( mesh )
    %         [ rx, ry, rz, sx, sy, sz, tx, ty, tz, J ] ...
    %             = mesh.cell.assembleJacobianMatrix( mesh.x, mesh.y, mesh.z );
    %     end
        
    %     function [ nx, ny, nz, Js ] = assembleNormalVector( mesh )
    %         [ nx, ny, nz, Js ] = mesh.cell.assembleNormalVector( mesh.x, mesh.y, mesh.z );
    %     end
        
    %     function LIFT = assembleLiftMatrix( cell )
    %         Mes = zeros(cell.Np, cell.TNfp);
    %         sk = 1;
    %         for f = 1:cell.Nface
    %             fcell = getStdCell(cell.N, cell.faceType(f));
    %             row = cell.Fmask(:, f);
    %             row = row(row ~= 0);
    %             for n = 1:fcell.Np
    %                 Mes(row, sk) = fcell.M(:, n);
    %                 sk = sk + 1;
    %             end
    %         end
    %         LIFT = cell.invM * Mes;
    %     end
    % end
    
end

