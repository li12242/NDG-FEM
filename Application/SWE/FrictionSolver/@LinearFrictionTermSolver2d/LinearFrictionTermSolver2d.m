classdef LinearFrictionTermSolver2d < AbstractFrictionTermSolver
    properties
        %> bottom friction coefficient
        r
    end
    
    methods
        function obj = LinearFrictionTermSolver2d(t)%���캯��
            obj.r = t;
        end
        
        function evaluateFrictionTermRHS( obj, physClass, fphys )

            a = obj.r;

            for m = 1:physClass.Nmesh
                
                mesh = physClass.meshUnion(m);
                ind = (mesh.EToR == int8(NdgRegionType.Wet));
                
                % frhs = frhs - rhu   
                physClass.frhs{m}(:,ind,2) = physClass.frhs{m}(:,ind,2)...
                    - a*(fphys{m}(:,ind,2));
                
                % frhs = frhs - rhv 
                physClass.frhs{m}(:,ind,3) = physClass.frhs{m}(:,ind,3)...
                    - a*(fphys{m}(:,ind,3));
                
            end
        end
    end
    
end
