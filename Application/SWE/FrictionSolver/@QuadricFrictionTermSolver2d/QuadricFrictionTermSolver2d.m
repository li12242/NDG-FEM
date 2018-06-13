classdef QuadricFrictionTermSolver2d < AbstractFrictionTermSolver
    
    properties
        %> bottom friction coefficient
        r
    end
    
    methods
        function obj = QuadricFrictionTermSolver2d(t)
            obj.r = t;
        end
        
        function evaluateFrictionTermRHS( obj, phys, fphys )
            
            for m = 1:phys.Nmesh
                
                mesh = phys.meshUnion(m);
                ind = (mesh.status == int8(enumSWERegion.Wet));
                
                s = sqrt( fphys{m}(:,ind,2).^2 + fphys{m}(:,ind,1).^2 )./fphys{m}(:,ind,1).^2;
                
                % frhs = frhs - rhu
                phys.frhs{m}(:,ind,2) = phys.frhs{m}(:,ind,2)...
                    - obj.r*(fphys{m}(:,ind,2) .* s);
                
                % frhs = frhs - rhv
                phys.frhs{m}(:,ind,3) = phys.frhs{m}(:,ind,3)...
                    - obj.r*(fphys{m}(:,ind,3) .* s);
                
            end
        end
    end
end