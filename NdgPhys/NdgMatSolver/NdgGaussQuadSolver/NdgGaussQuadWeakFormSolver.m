classdef NdgGaussQuadWeakFormSolver < NdgGaussQuadStrongFormSolver
    
    methods
        function obj = NdgGaussQuadWeakFormSolver( phys )
            obj = obj@NdgGaussQuadStrongFormSolver( phys );

            for m = 1:phys.Nmesh
                mesh = phys.meshUnion( m );
                [ obj.Vq{m} ] = mesh.cell.Vq;
                obj.Dr{m} = obj.Dr{m}';
                obj.Ds{m} = obj.Ds{m}';
                obj.Dt{m} = obj.Dt{m}';
            end
        end% func
    end
    
end

