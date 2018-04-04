classdef SWEWD1d < SWEPreBlanaced1d
    
    methods
        function [ E ] = matEvaluateFlux( obj, mesh, fphys )
            [ E ] = mxEvaluateFlux1d( obj.hmin, obj.gra, mesh.EToR, fphys );
        end
    end
    
    methods( Access = protected )
        
        function matEvaluateTopographySourceTerm( obj, fphys )
            for m = 1:obj.Nmesh
                mesh = obj.meshUnion(m);
                obj.frhs{m} = obj.frhs{m} + mxEvaluateSourceTopography1d...
                    ( obj.gra, mesh.EToR, fphys{m}, obj.zGrad{m} );
            end
        end
        
        function [ fphys ] = matEvaluatePostFunc(obj, fphys)
            obj.matUpdateWetDryState( fphys );
        end% func
        
        function fphys = matEvaluateLimiter( obj, fphys )
            fphys = matEvaluateLimiter@SWEAbstract1d( obj, fphys );
            obj.matUpdateWetDryState( fphys )
            for m = 1:obj.Nmesh
                mesh = obj.meshUnion(m);
                
%                 ind1 = temp(1, :) >= 0;
%                 ind2 = ~ind1;
%                 pos = sum( temp(1, ind1) );
%                 neg = sum( temp(1, ind2) );
%                 fprintf('pos = %e, neg = %e, tol = %12.10e\n', pos, neg, pos-neg)
%                 temp(1, ind1) = temp(1, ind1) + temp(1, ind1)*neg/pos;
%                 temp(:, ind2) = 0;
                
                vol = mesh.GetMeshIntegralValue( fphys{m}(:, :, 1) );
                ind1 = (vol >= 0);
                pos = sum( vol( ind1 ) );
                neg = sum( vol( ~ind1 ) );
                fphys{m}(:, ind1, 1) = fphys{m}(:, ind1, 1)*( 1 + neg/pos );
                fphys{m}(:, ~ind1, 1) = 0.0;
                
                flg = (mesh.EToR == int8( NdgRegionType.PartialWetDamBreak ) ) ...
                    | (mesh.EToR == int8( NdgRegionType.PartialWetFlood ));
                temp = mesh.cell.V \ fphys{m}(:, flg, 1);
                temp(3:end, :) = 0;
                fphys{m}(:, flg, 1) = mesh.cell.V * temp;
                
%                 for fld = 2:obj.Nvar % reconstruct the WD cell values
%                     tempWD = mesh.cell.V \ fphys{m}(:, flg, fld);
%                     ind = tempWD(1, :) < 0;
%                     tempWD(1, ind) = 0;
%                     tempWD(3:end, :) = 0;
%                     fphys{m}(:, flg, fld) = mesh.cell.V * tempWD;
%                 end
            end    
        end% func
                
        function matUpdateWetDryState(obj, fphys)
            for m = 1:obj.Nmesh
                
                mesh = obj.meshUnion(m);
                % wet part
                wetflg = all( fphys{m}(:,:,1) > obj.hmin );
                mesh.EToR( wetflg ) = int8( NdgRegionType.Wet );
                % dry part
                dryflag = all( fphys{m}(:,:,1) < obj.hmin );
                mesh.EToR( dryflag ) = int8( NdgRegionType.Dry );
                % partial wet part
                pwlg = ( ~wetflg ) & ( ~dryflag ); %Nwd = sum( flg );
                % dram-break type or flood type
                dblg = max( fphys{m}(:,:,1) + fphys{m}(:,:,3) ) ...
                    > max( fphys{m}(:,:,3) );
                pdblg = dblg & pwlg;
                pflg = (~dblg) & pwlg;
                mesh.EToR( pdblg ) = int8( NdgRegionType.PartialWetDamBreak );
                mesh.EToR( pflg ) = int8( NdgRegionType.PartialWetFlood );
                %mesh.EToR( flg ) = int8( NdgRegionType.PartialWet );
            end
        end% func
    end
    
end

