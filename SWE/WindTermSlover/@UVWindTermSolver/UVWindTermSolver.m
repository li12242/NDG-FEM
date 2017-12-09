classdef UVWindTermSolver < AbstractWindTermSolver
    
    properties
        cd
        rou
        rouair
    end
    
    methods
        function obj = StressWindTermSolver(o, p, q)
            obj.cd = o;
            obj.rou = p;
            obj.rouair = q;
        end
        
        function evaluateWindTermRHS( obj, physClass, fphys )
            %> wind stress coefficient
            b = obj.cd;
            %> dencity of ocean(kg/m^3)
            c = obj.rou;
            %> dencity of air(kg/m^3)
            d = obj.rouair;
            
            for m = 1:physClass.Nmesh
                
                mesh = physClass.meshUnion(m);
                ind = (mesh.EToR == int8(NdgRegionType.Wet));

                
                w10 = sqrt(fphys{m}(:,:,7).*fphys{m}(:,:,7)+fphys{m}(:,:,8).*fphys{m}(:,:,8));
                
                % frhs = frhs + cd*rouair*w10*windx/rouwater
                physClass.frhs{m}(:,ind,2) = physClass.frhs{m}(:,ind,2)...
                    + (b*d*w10(:,ind).*(fphys{m}(:,ind,7)))/c;
                
                % frhs = frhs + cd*rouair*w10*windx/rouwater
                physClass.frhs{m}(:,ind,3) = physClass.frhs{m}(:,ind,3)...
                    + (b*d*w10(:,ind).*(fphys{m}(:,ind,8)))/c;
                
            end
        end
    end
    
end

