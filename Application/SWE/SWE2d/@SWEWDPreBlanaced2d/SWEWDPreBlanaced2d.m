classdef SWEWDPreBlanaced2d < SWEPreBlanaced2d
    
    methods( Access = protected )
        
        function matUpdateWDWetDryState(obj, fphys)
            for n = 1:obj.Nmesh
                mesh = obj.meshUnion(n);
                mesh.EToR = mxUpdateWDWetDryState( obj.hmin, fphys{n} );
            end
        end% func
        
        function [ fphys ] = matEvaluatePostFunc(obj, fphys)
            for m = 1:obj.Nmesh
                mesh = obj.meshUnion(m);
                flg = ( mesh.EToR == int8( enumSWERegion.PartialWet ) );
                
                tempWD = mesh.cell.V \ fphys{m}(:, flg, 1);
                ind = tempWD(1, :) < 0;
                tempWD(1, ind) = 0;
                tempWD(4:end, :) = 0;
                fphys{m}(:, flg, 1) = mesh.cell.V * tempWD;
                
                for fld = 2:obj.Nvar % reconstruct the WD cell values
                    tempWD = mesh.cell.V \ fphys{m}(:, flg, fld);
                    tempWD(4:end, :) = 0;
                    fphys{m}(:, flg, fld) = mesh.cell.V * tempWD;
                end
            end
            
            obj.matUpdateWetDryState( fphys );
        end% func
        
%         function fphys = matEvaluateLimiter( obj, fphys )
%             obj.matUpdateWDWetDryState( fphys )
%             fphys = matEvaluateLimiter@SWEPreBlanaced2d( obj, fphys );
%         end% func
    end
    
end

