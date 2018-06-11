classdef BetaApproCoriolisTermSolver < AbstractCoriolisTermSolver
    properties
        f0
        beta
    end
    
    methods
        function obj = BetaApproCoriolisTermSolver(m, n)
            obj.f0 = m;
            obj.beta = n;
        end
        
        function evaluateCoriolisTermRHS( obj, physClass, fphys )

            a = obj.f0;
            b = obj.beta;

            for m = 1:physClass.Nmesh 
                
                mesh = physClass.meshUnion(m);
                ind = (mesh.EToR == int8(NdgRegionType.Wet));

    
                Np = physClass.meshUnion(m).cell.Np;
                k = physClass.meshUnion(m).K;
                s = ones(Np,k);
                q = a*s;%f0
                
                % frhs = frhs + (f+by)hv
                physClass.frhs{m}(:,ind,2) = physClass.frhs{m}(:,ind,2)...
                    + (q(:,ind)+b*mesh.y(:,ind)).*(fphys{m}(:,ind,3));
                
                % frhs = frhs - (f+by)hu
                physClass.frhs{m}(:,ind,3) = physClass.frhs{m}(:,ind,3)...
                    - (q(:,ind)+b*mesh.y(:,ind)).*(fphys{m}(:,ind,2));
                
            end
        end
    end
    
end

