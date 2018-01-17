classdef SWEWDMesh1d < SWEConventional1d
    
    methods( Hidden )
        function initPhysFromOptions( obj, mesh )
            EToB = [1, mesh.Nv; double( mesh.eidtype([1, end]) )];
            mesh2 = NdgMesh1d( ...
                StdLine(1), mesh.Nv, mesh.vx, ...
                mesh.K, mesh.EToV, mesh.EToR, EToB );
            mesh = [mesh, mesh2];
            initPhysFromOptions@SWEAbstract1d( obj, mesh );
            %obj.advectionSolver = NdgAdvWDSolver1d( obj );
        end
        
        function [ fM, fP ] = matEvaluateSurfaceValue( obj, mesh, fphys, fext )
            [ fM, fP ] = mxEvaluateSurfaceValue1d( obj.hmin, obj.gra, ...
                mesh.EToR, mesh.eidM, mesh.eidP, ...
                mesh.nx, mesh.eidtype, fphys, fext );
        end
    end
    
    methods( Access = protected )
        
        function fphys = matUpdateExternalField( obj, tloc, fphys )
            % update the external value of low order field to high order
            % field
            mesh = obj.meshUnion(1);            
            flg = (mesh.EToR == int8( NdgRegionType.PartialWet ) );
            meshWd = obj.meshUnion(2);
%             Nwd = sum( flg );
            for fld = 1:obj.Nvar % reconstruct the WD cell values
                temp = mesh.cell.V \ fphys{1}(:, flg, fld);
                tempWD = temp(1:2, :);
                fphys{2}(:, flg, fld) = meshWd.cell.V * tempWD;
                
                obj.fext{1}([1,end], flg, fld) = fphys{2}([1,end], flg, fld);
                obj.fext{2}([1,end], flg, fld) = fphys{1}([1,end], flg, fld);
                
            end
        end
        
        function fphys = matEvaluateLimiter( obj, fphys )
            mesh = obj.meshUnion(1);
            meshWd = obj.meshUnion(2);
            
            flg = (mesh.EToR == int8( NdgRegionType.PartialWet ) );
            Nwd = sum( flg );
            for fld = 1:obj.Nvar % reconstruct the WD cell values
                tempWD = meshWd.cell.V \ fphys{2}(:, flg, fld);
                temp = zeros( mesh.cell.Np, Nwd );
                temp(1:2, :) = tempWD;
                fphys{1}(:, flg, fld) = mesh.cell.V * temp;
                
                temp = mesh.cell.V \ fphys{1}(:, ~flg, fld);
                tempWD = temp(1:2, :);
                fphys{2}(:, ~flg, fld) = meshWd.cell.V * tempWD;
            end
            
            fphys = matEvaluateLimiter@SWEAbstract1d( obj, fphys );
        end
        
        function [ fphys ] = matEvaluatePostFunc(obj, fphys)
            obj.matUpdateWetDryState( fphys );
        end
        
        function matUpdateWetDryState(obj, fphys)
            mesh = obj.meshUnion(1);
            meshWd = obj.meshUnion(2);            
            
            % wet part
            wetflg = all( fphys{1}(:,:,1) > obj.hmin );
            mesh.EToR( wetflg ) = int8( NdgRegionType.Wet );
            meshWd.EToR( wetflg ) = int8( NdgRegionType.Wet );
            % dry part
            dryflag = all( fphys{1}(:,:,1) < obj.hmin );
            mesh.EToR( dryflag ) = int8( NdgRegionType.Dry );
            meshWd.EToR( dryflag ) = int8( NdgRegionType.Dry );
            
            % partial wet part
            flg = ( ~wetflg ) & ( ~dryflag ); %Nwd = sum( flg );
            mesh.EToR( flg ) = int8( NdgRegionType.PartialWet );
            meshWd.EToR( flg ) = int8( NdgRegionType.PartialWet );

        end% func
    end
    
end

