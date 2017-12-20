classdef NdgQuadFreeWeakFormSlver < NdgQuadFreeStrongFormSolver
    %NDGQUADFREEWEAKFORMSLVER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        M
    end
    
    methods
        function obj = NdgQuadFreeWeakFormSlver( phys )
            obj = obj@NdgQuadFreeStrongFormSolver( phys );
            
            for m = 1:phys.Nmesh
                mesh = phys.meshUnion( m );
                [ obj.Dr{m}, obj.Ds{m}, obj.Dt{m} ] = obj.assembleDerivativeMatrix( mesh.cell );
                [ obj.M{m} ] = mesh.cell.M;
            end
        end
    end
    
    methods( Static )
        function [Dr, Ds, Dt] = assembleDerivativeMatrix( cell )
            Vr = zeros(cell.Np, cell.Np);
            Vs = zeros(cell.Np, cell.Np);
            Vt = zeros(cell.Np, cell.Np);
            for n = 1:cell.Np
                [Vr(:, n), Vs(:, n), Vt(:, n)] = cell.orthogonal_derivative_func...
                    (n, cell.r, cell.s, cell.t);
            end
            Dr = (cell.V * Vr');
            Ds = (cell.V * Vs');
            Dt = (cell.V * Vt');
        end
    end
    
end

