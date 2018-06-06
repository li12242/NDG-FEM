classdef ManningFrictionSolver1d < AbstractFrictionTermSolver
    %MANNINGFRICTIONSOLVER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        n
    end
    
    methods
        function obj = ManningFrictionSolver2d( n )
            obj.n = n;
        end% func
        
        function evaluateFrictionTermRHS( obj, physObj, fphys )
            for m = 1:physObj.Nmesh
                
                mesh = physObj.meshUnion(m);
                ind = (mesh.EToR == int8(NdgRegionType.Wet));
                
                qn  = abs( fphys{m}(:,ind,2) );
                % frhs = frhs - rhu   
                physObj.frhs{m}(:,ind,2) = physObj.frhs{m}(:,ind,2)...
                    - obj.n .* physObj.gra * fphys{m}(:,ind,2) .* qn ...
                    ./( fphys{m}(:,ind,1).^(7/3) ) ;
            end
        end% func
    end
    
end

