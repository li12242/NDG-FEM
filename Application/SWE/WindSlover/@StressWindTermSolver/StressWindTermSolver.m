classdef StressWindTermSolver < AbstractWindTermSolver
    
    properties
        rou
    end
    
    methods
        function obj = StressWindTermSolver(q)
            obj.rou = q;
        end
        
        function evaluateWindTermRHS( obj, physClass, fphys )

            %> dencity of ocean(kg/m^3)
            c = obj.rou;

            for m = 1:physClass.Nmesh
                
                mesh = physClass.meshUnion(m);
                ind = (mesh.EToR == int8(NdgRegionType.Wet));

                
                % frhs = frhs + tx/rouwater
                physClass.frhs{m}(:,ind,2) = physClass.frhs{m}(:,ind,2)...
                    + fphys{m}(:,ind,7)/c;
                
                % frhs = frhs + ty/rouwater
                physClass.frhs{m}(:,ind,3) = physClass.frhs{m}(:,ind,3)...
                    + fphys{m}(:,ind,8)/c;
                
            end
        end
    end
    
end

