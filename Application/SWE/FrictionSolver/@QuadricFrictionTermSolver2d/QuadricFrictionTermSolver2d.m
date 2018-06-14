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
            
            a = obj.r;
            g = phys.gra;
            
            for m = 1:phys.Nmesh
                
                mesh = phys.meshUnion(m);
                ind = (mesh.status == int8(enumSWERegion.Wet));
                
                s = sqrt( fphys{m}(:,ind,2).^2 + fphys{m}(:,ind,3).^2 )./fphys{m}(:,ind,1).^2;
                
                % frhs = frhs - rhu
                phys.frhs{m}(:,ind,2) = phys.frhs{m}(:,ind,2)...
                    - g*a*a*(fphys{m}(:,ind,2) .* s)./(fphys{m}(:,ind,1).^(1/3));   
                
                % frhs = frhs - rhv
                phys.frhs{m}(:,ind,3) = phys.frhs{m}(:,ind,3)...
                    - g*a*a*(fphys{m}(:,ind,3) .* s)./(fphys{m}(:,ind,1).^(1/3));
                
            end
        end
    end
end