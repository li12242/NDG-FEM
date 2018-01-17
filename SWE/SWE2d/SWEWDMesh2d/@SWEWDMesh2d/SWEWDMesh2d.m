classdef SWEWDMesh2d < SWEPreBlanaced2d
    %SWEWDMESH2D Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods( Hidden )
        function initPhysFromOptions( obj, mesh )
            EToB = [1, mesh.Nv; double( mesh.eidtype([1, end]) )];
            mesh2 = NdgMesh1d( ...
                StdLine(1), mesh.Nv, mesh.vx, ...
                mesh.K, mesh.EToV, mesh.EToR, EToB );
            mesh = [mesh, mesh2];
            initPhysFromOptions@SWEPreBlanaced1d( obj, mesh );
            %obj.advectionSolver = NdgAdvWDSolver1d( obj );
        end
    end
    
    methods( Access = protected )
        
        function matUpdateExternalField( obj, tloc, fphys )
            % update the external value of low order field to high order
            % field
            mesh = obj.meshUnion(1);
            meshWd = obj.meshUnion(2);
            flg = (mesh.EToR == int8( NdgRegionType.PartialWet ) );
            Nwd = sum( flg );
            for fld = 1:obj.Nfield % reconstruct the WD cell values
                tempWD = meshWd.cell.V \ fphys{2}(:, flg, fld);
                temp = zeros( mesh.cell.Np, Nwd );
                temp(1:3, :) = tempWD;
                obj.fext{1}(:, flg, fld) = mesh.cell.V * temp;
                
                temp = mesh.cell.V \ fphys{1}(:, flg, fld);
                tempWD = temp(1:3, :);
                obj.fext{2}(:, flg, fld) = meshWd.cell.V * tempWD;
            end
        end
        
%         function [ fphys ] = matEvaluatePostFunc(obj, fphys)
%             obj.matUpdateWetDryState( fphys );
%         end
        
        function matUpdateWetDryState(obj, fphys)
            mesh = obj.meshUnion(1);
            meshWd = obj.meshUnion(2);
            
            flg = (mesh.EToR == int8( NdgRegionType.PartialWet ) );
            Nwd = sum( flg );
            for fld = 1:obj.Nfield % reconstruct the WD cell values
                tempWD = meshWd.cell.V \ fphys{2}(:, flg, fld);
                ind = (tempWD(1, :) < 0);
                tempWD(1, ind) = 0;
                temp = zeros( mesh.cell.Np, Nwd );
                temp(1:3, :) = tempWD;
                fphys{1}(:, flg, fld) = mesh.cell.V * temp;
                
                temp = mesh.cell.V \ fphys{1}(:, ~flg, fld);
                tempWD = temp(1:3, :);
                fphys{2}(:, ~flg, fld) = meshWd.cell.V * tempWD;
%                 tempWD = meshWd.cell.V \ resQ{2}(:, flg, fld);
%                 temp(1:2, :) = tempWD;
%                 resQ{1}(:, flg, fld) = mesh.cell.V * temp;
            end
            
            % wet part
            wetflg = all( fphys{1}(:,:,1) > obj.hmin );
            mesh.EToR( wetflg ) = int8( NdgRegionType.Wet );
            % dry part
            dryflag = all( fphys{1}(:,:,1) < obj.hmin );
            mesh.EToR( dryflag ) = int8( NdgRegionType.Dry );
            fphys{1}(:, dryflag, 2) = 0;
            fphys{2}(:, dryflag, 2) = 0;
            
            % partial wet part
            flg = ( ~wetflg ) & ( ~dryflag );
            mesh.EToR( flg ) = int8( NdgRegionType.PartialWet );
            for fld = 1:obj.Nfield % reconstruct the WD cell values
                temp = mesh.cell.V \ fphys{1}(:, flg, fld);
                tempWD = temp(1:3, :);
                fphys{2}(:, flg, fld) = meshWd.cell.V * tempWD;
                
%                 temp = mesh.V \ resQ{1}(:, flg, fld);
%                 tempWD = temp(1:2, :);
%                 resQ{2}(:, flg, fld) = meshWd.cell.V * tempWD;
            end
        end% func
    end
end

