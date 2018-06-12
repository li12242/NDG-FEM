classdef SWEWDPreBlanaced2d < SWEPreBlanaced2d
    
    methods( Hidden )
        function [ E, G ] = matEvaluateFlux( obj, mesh, fphys )
            [ E, G ] = mxEvaluateFlux2d( obj.hmin, obj.gra, mesh.status, fphys );
        end
    end
    
    methods( Access = protected )
        
%         function [ fphys ] = matEvaluatePostFunc(obj, fphys)
%             for m = 1:obj.Nmesh
%                 mesh = obj.meshUnion(m);
%                 flg = ( mesh.status == ...
%                     int8( enumSWERegion.PartialWetDamBreak ) ) ...
%                     | ( mesh.status == ...
%                     int8( enumSWERegion.PartialWetFlood ) );
%                 
%                 tempWD = mesh.cell.V \ fphys{m}(:, flg, 1);
%                 tempWD(1, :) = max( 0, tempWD(1, :) );
%                 tempWD(4:end, :) = 0;
%                 fphys{m}(:, flg, 1) = mesh.cell.V * tempWD;
%                 
%                 tempWD = mesh.cell.V \ fphys{m}(:, flg, 2);
%                 tempWD(4:end, :) = 0;
%                 fphys{m}(:, flg, 2) = mesh.cell.V * tempWD;
%                 
%                 tempWD = mesh.cell.V \ fphys{m}(:, flg, 3);
%                 tempWD(4:end, :) = 0;
%                 fphys{m}(:, flg, 3) = mesh.cell.V * tempWD;
%             end
%             obj.matUpdateWetDryState( fphys );
%         end% func
        
        function matUpdateWetDryState(obj, fphys)
            for n = 1:obj.Nmesh
                mesh = obj.meshUnion(n);
                mesh.status = mxUpdateWDWetDryState( obj.hmin, fphys{n} );
            end
        end% func
        
        function matEvaluateTopographySourceTerm( obj, fphys )
            for m = 1:obj.Nmesh
                mesh = obj.meshUnion(m);
                obj.frhs{m} = obj.frhs{m} + mxEvaluateSourceTopography2d...
                    ( obj.gra, mesh.status, fphys{m}, obj.zGrad{m} );
            end
        end
    end
    
end

