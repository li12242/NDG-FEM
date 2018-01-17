classdef SWESubcellPreBalance1d < SWEConventional1d
    
    properties( SetAccess = protected )
        %> mesh contains the wet/dry elements
        meshWD
    end
    
    methods
        function initPhysFromOptions( obj, mesh )
            initPhysFromOptions@SWEConventional1d( obj, mesh );
            obj.meshWD = NdgWDMesh1d( mesh ); % init subcell mesh
            obj.advectionSolver = NdgAdvSubWDSolver1d( obj );
            obj.Nmesh = 1;
            mesh = obj.meshWD;
            obj.fphys{2} = zeros( mesh.cell.Np, mesh.K, obj.Nfield );
            obj.frhs{2} = zeros( mesh.cell.Np, mesh.K, obj.Nvar );
            obj.fext{2} = zeros( mesh.cell.Np, mesh.K, obj.Nfield );
        end
        
        function drawSubMesh( obj, varargin )
            if nargin > 1
                fphys = varargin{1}{2};
                obj.draw( varargin{1} );
            else
                fphys = obj.fphys{2};
                obj.draw( obj.fphys );
            end
            
            mesh = obj.meshWD;
            K = mesh.K;
            if K > 0
                if numel( obj.draw_handle ) == 4
                    obj.draw_handle{2,1}.delete;
                    obj.draw_handle{2,2}.delete;
                end
                list = 1:mesh.cell.Np:(K*mesh.cell.Np);
                g = graph();
                for n = 1:mesh.cell.Np-1
                    g = addedge(g, list+n-1, list+n);
                end
                subplot(2,1,1); hold on;
                x = mesh.x(:, 1:K);
                temp = fphys(:, 1:K, 1) + fphys(:, 1:K, 3);
                obj.draw_handle{2,1} = plot(g, ...
                    'XData', x(:), 'YData', temp(:), ...
                    'LineWidth', 1, ...
                    'Marker', 'o', ...
                    'NodeColor','g', ...
                    'EdgeColor', 'g', ...
                    'MarkerSize', 2, ...
                    'NodeLabel', {});
                subplot(2,1,2); hold on;
                temp = fphys(:, 1:K, 2);
                obj.draw_handle{2,2} = plot(g, ...
                    'XData', x(:), 'YData', temp(:), ...
                    'LineWidth', 1, ...
                    'Marker', 'o', ...
                    'NodeColor','g', ...
                    'EdgeColor', 'g', ...
                    'MarkerSize', 2, ...
                    'NodeLabel', {});
                
            end
        end% function
    end
    
    methods( Access = protected )
        matEvaluateRK45( obj );
        
        function matEvaluateTopographySourceTerm( obj, fphys )
        end
        
        function matUpdateExternalField( obj, time, fphys )
            matUpdateExternalField@NdgPhysMat( obj, time, fphys );
            for k = 1:(obj.meshWD.K/2)
                k1 = 2*k - 1;
                k2 = 2*k;
                cid = obj.meshWD.coarseCellId(k1);
                eid = obj.meshUnion(1).EToE(:, cid);
                obj.fext{1}(end, eid(1), :) = fphys{2}(1, k1, :);
                obj.fext{1}(1, eid(2), :) = fphys{2}(end, k2, :);
                
                obj.fext{2}([1,end], k1, :) = fphys{1}([1,end], eid(1), :);
                obj.fext{2}([1,end], k2, :) = fphys{1}([1,end], eid(2), :);
            end
        end% function   
        %>
        function [ fphys ] = matEvaluatePostFunc(obj, fphys)
            fphys = obj.recombineWDMesh( fphys );
            if obj.isAnyNan( fphys )
                keyboard;
            end
            obj.matUpdateWetDryState( fphys );
            fphys = obj.resetWDMesh( fphys );
            if obj.isAnyNan( fphys )
                keyboard;
            end
        end
        
        function flag = isAnyNan( obj, fphys )
            fld = 1;
            flag = any( any( isnan(fphys{1}(:,:,fld)) ) ) | ...
                any( any( isnan(fphys{2}(:,:,fld)) ) );
            fld = 2;
            flag = flag || any( any( isnan(fphys{1}(:,:,fld)) ) ) || ...
                any( any( isnan(fphys{2}(:,:,fld)) ) );
        end
        
        function matUpdateWetDryState(obj, fphys)
            for m = 1:obj.Nmesh
                wetflag = all( fphys{m}(:,:,1) > obj.hmin );
                obj.meshUnion(m).EToR( ~wetflag ) = int8( NdgRegionType.Dry );
                obj.meshUnion(m).EToR(  wetflag ) = int8( NdgRegionType.Wet );
                pwetflag = any( fphys{m}(:,:,1) > obj.hmin ) & ~wetflag;
                obj.meshUnion(m).EToR(  pwetflag ) = int8( NdgRegionType.PartialWet );
            end
        end
    end
    
    methods( Access = private )
        
        function fphys = recombineWDMesh( obj, fphys )
            mesh = obj.meshUnion(1);
            K = obj.meshWD.K;
            hc = obj.meshWD.GetMeshIntegralValue( ...
                obj.meshWD.J(:,1:K), obj.meshWD.LAV(1:K), fphys{2}(:, 1:K, 1) );
            qc = obj.meshWD.GetMeshIntegralValue( ...
                obj.meshWD.J(:,1:K), obj.meshWD.LAV(1:K), fphys{2}(:, 1:K, 2) );
            for k = 1:(K/2)
                k1 = 2*k - 1;
                k2 = 2*k;
                c1 = obj.meshWD.coarseCellId(k1);
                dx1 = abs( diff(obj.meshWD.x( [1,end], k1 )) );
                dx2 = abs( diff(obj.meshWD.x( [1,end], k2 )) );
                dxc = obj.meshWD.xc(k2) - obj.meshWD.xc(k1);
                slope = (hc(k2) - hc(k1))./dxc;
                hm = max(0, (hc(k1) * dx1 + hc(k2) * dx2)/ mesh.LAV(c1));
                ratio = min( 2*abs(hm) ./ mesh.LAV(c1)./abs(slope), 1 );
%                 fphys{1}(:, c1, 1) = hm;
                dx = mesh.x(:, c1) - mesh.xc(c1);
                fphys{1}(:, c1, 1) = bsxfun(@times, slope.*ratio, dx ) + hm;
                
                slope = (qc(k2) - qc(k1))./dxc;
                qm = (qc(k1) * dx1 + qc(k2) * dx2)/ mesh.LAV(c1);
                ratio = min( 2*abs(qm) ./ mesh.LAV(c1)./abs(slope), 1 );
%                 fphys{1}(:, c1, 2) = qm;
                fphys{1}(:, c1, 2) = bsxfun(@times, slope.*ratio, dx ) + qm;
                % reset coarse edge type
                eid = mesh.EToE(:, c1);
                fid = mesh.EToF(:, c1);
                nid = sub2ind([2, mesh.K], fid, eid);
                mesh.eidtype( nid ) = NdgEdgeType.Inner;
            end
        end% function
        
        function fphys = resetWDMesh( obj, fphys )
            mesh = obj.meshUnion(1);
            meshWDId = 2;
            obj.meshWD.reinitSubcellMesh( meshWDId );
            ids = find( mesh.EToR == int8( NdgRegionType.PartialWet ) );
            Nwd = numel( ids ); % number of partial wet cell
            hc = mesh.GetMeshAverageValue( fphys{1}(:,:,1) );
            qc = mesh.GetMeshAverageValue( fphys{1}(:,:,2) );
            for n = 1:Nwd
                cid = ids(n);
                b = fphys{1}(:,cid,3);
                h = fphys{1}(:,cid,1);
                eta = h + b;
                % reconstruct water depth and flow flux
                if any( eta > b ) % dam break type
                    [ hpmax ] = max( h( [1,end] ) );
                    if hpmax - 2*hc(cid) < 1e-5
                        [ fphys ] = obj.reconstructDambreakStates( ...
                            mesh, fphys, cid);
                        continue;
                    end
                    [xs, fphys, cellType] = obj.reconstructDamBreakWDState( ...
                        mesh, fphys, hc(cid), qc(cid), cid, n);
                else% flood type
                    [xs, fphys, cellType] = obj.reconstructFloodWDState(...
                        mesh, fphys, hc(cid), qc(cid), cid, n);
                end
                obj.meshWD.setSubcell([mesh.x( [1,end], cid ); xs], mesh.eidP(:, cid), ...
                    cellType, ...
                    [NdgEdgeType.RefineEdge, NdgEdgeType.RefineEdge], meshWDId, cid);
                % set the adjacent element edge type
                eid = mesh.EToE(:, cid);
                fid = mesh.EToF(:, cid);
                nid = sub2ind([2, mesh.K], fid, eid);
                mesh.eidtype( nid ) = NdgEdgeType.RefineEdge;
            end
        end% function
        
        function [xs, fphys, cellType] = reconstructFloodWDState(...
                obj, mesh, fphys, hc, qc, cellId, n )
            cell = obj.meshWD.cell;
            b = fphys{1}(:,cellId,3);
            xb = mesh.x( [1,end], cellId );
            dx = mesh.LAV( cellId );
            db = abs( b(end) - b(1) );
            [ ~, nid ] = min( b([1,end]) );
            subid = subcellId*2 - 1; % start subcell index
            switch nid
                case 1 % left is wet
                    xs = xb(1) + sqrt( 2*hc*dx.^2/db );
                    hmax = db / dx * ( xs - xb(1) );
                    fphys{2}(:, subid  , 1) = cell.project_vert2node( [hmax, 0]' );
                    fphys{2}(:, subid+1, 1) = 0;
                    cellType = [NdgRegionType.Wet, NdgRegionType.Dry];
                case 2 % right is wet
                    xs = xb(2) - sqrt( 2*hc*dx.^2/db );
                    hmax = db / dx * ( xb(2) - xs );
                    fphys{2}(:, subid  , 1) = 0;
                    fphys{2}(:, subid+1, 1) = cell.project_vert2node( [0, hmax]' );
                    cellType = [NdgRegionType.Dry, NdgRegionType.Wet];
            end
            fphys{2}(:, subid  , 2) = 0;
            fphys{2}(:, subid+1, 2) = 0;
            [ fphys{2}(:, subid, 3), fphys{2}(:, subid+1, 3) ] ...
                = obj.reconstructSubBottom( mesh, fphys, xs, cellId );
        end
        
        function [xs, fphys,cellType] = reconstructDamBreakWDState( ...
                obj, mesh, fphys, hc, qc, cellId, n )
            
            cell = obj.meshWD.cell;
            [ hpmax, nid ] = max( fphys{1}([1,end],cellId,1) );
            xb = mesh.x( [1,end], cellId );
            lamda = 2*hc/hpmax;
            subid = n*2 - 1; % start subcell index
            switch nid
                case 1
                    xs = lamda * xb( 2 ) + ( 1 - lamda ) * xb( 1 );
                    fphys{2}(:, subid, 1) = cell.project_vert2node( [hpmax, 0]' );
                    fphys{2}(:, subid, 2) = cell.project_vert2node...
                        ( [2*qc, 0]' ); %( [2*qc*( abs(xb(2) - xb(1)) )/abs( xb(1) - xs ), 0]' );
                    fphys{2}(:, subid+1, 1:2) = 0;
                    cellType = [NdgRegionType.Wet, NdgRegionType.Dry];
                case 2
                    xs = lamda * xb( 1 ) + ( 1 - lamda ) * xb( 2 );
                    fphys{2}(:, subid+1, 1) = cell.project_vert2node( [0, hpmax]' );
                    fphys{2}(:, subid+1, 2) = cell.project_vert2node...
                        ( [0, 2*qc]' ); %( [0, 2*qc*( dx )/abs( xb(1) - xs )]' );
                    fphys{2}(:, subid, 1:2) = 0;
                    cellType = [NdgRegionType.Dry, NdgRegionType.Wet];
            end
            fphys{2}(:, subid  , 2) = 0;
            fphys{2}(:, subid+1, 2) = 0;
            [ fphys{2}(:, subid, 3), fphys{2}(:, subid+1, 3) ] ...
                = obj.reconstructSubBottom( mesh, fphys, xs, cellId );
        end
        
        function [ b1, b2 ] = reconstructSubBottom( obj, mesh, fphys, xs, cid )
            xb = mesh.x( [1,end], cid );
            bot = fphys{1}( [1,end], cid, 3 );
            bs = bot(1) + ( bot(2) - bot(1) ) * ( xs - xb(1) )./( xb(2) - xb(1) );
            b1 = obj.meshWD.cell.project_vert2node( [bot(1), bs]' );
            b2 = obj.meshWD.cell.project_vert2node( [bs, bot(2)]' );
        end
        
        %> Utilize linear function to reconstruct the dam break type transition cell.
        function [ fphys ] = reconstructDambreakStates( obj, mesh, fphys, cid )
            fm = mesh.cell.V\fphys{1}(:, cid, 1);
            temp = zeros( size(fm) );
            temp(1:2) = fm(1:2);
            fphys{1}(:, cid, 1) = mesh.cell.V*temp;
            
            qm = mesh.cell.V\fphys{1}(:, cid, 2);
            temp = zeros( size(fm) );
            temp(1:2) = qm(1:2);
            fphys{1}(:, cid, 2) = mesh.cell.V*temp;
            mesh.EToR( cid ) = NdgRegionType.Wet;
        end
        
%         function [ xs, fphys ] = reconstructWDState( obj, mesh, fphys, ...
%                 hc, qc, cellId, subcellId )
% 
%             cell = obj.meshUnion.cell;
%             h = fphys{1}(:,cellId,1);
%             b = fphys{1}(:,cellId,3);
%             eta = h + b;
%             xb = mesh.x( [1,end], cellId );
%             dx = mesh.LAV( cellId );
%             if max( eta ) > max( b ) % dambreak type
%                 %hb = h( cell.Fmask(:) );
%                 [ hpmax, nid ] = max( h( [1,end] ) );
%                 hpmax = max( hpmax, 2*hc );
% %                 if hpmax < 2*hc
% %                     warning([ 'The WD treatment in #', num2str(cellId),...
% %                         ' cell is incorrect.']);
% %                 end
%                 lamda = 2*hc/hpmax;
%                 subid = subcellId*2 - 1; % start subcell index
%                 switch nid
%                     case 1 % left cell is wet
%                         xs = lamda * xb( 2 ) + ( 1 - lamda ) * xb( 1 );
%                         fphys{2}(:, subid, 1) = cell.project_vert2node( [hpmax, 0]' );
%                         fphys{2}(:, subid, 2) = cell.project_vert2node...
%                             ( [2*qc, 0]' );
%                             %( [2*qc*( abs(xb(2) - xb(1)) )/abs( xb(1) - xs ), 0]' );
%                         fphys{2}(:, subid+1, 1:2) = 0;
%                     case 2 % right cell is wet
%                         xs = lamda * xb( 1 ) + ( 1 - lamda ) * xb( 2 );
%                         fphys{2}(:, subid+1, 1) = cell.project_vert2node( [0, hpmax]' );
%                         fphys{2}(:, subid+1, 2) = cell.project_vert2node...
%                             ( [0, 2*qc]' );
%                             %( [0, 2*qc*( dx )/abs( xb(1) - xs )]' );
%                         fphys{2}(:, subid, 1:2) = 0;
%                 end
% 
%             else % flood type
%                 db = abs( b(end) - b(1) );
%                 [ ~, nid ] = min( b([1,end]) );
%                 subid = subcellId*2 - 1; % start subcell index
%                 switch nid
%                     case 1 % left is wet
%                         xs = xb(1) + sqrt( 2*hc*dx.^2/db );
%                         hmax = db / dx * ( xs - xb(1) );
%                         fphys{2}(:, subid  , 1) = cell.project_vert2node( [hmax, 0]' );
%                         fphys{2}(:, subid+1, 1) = 0;
%                         
%                     case 2 % right is wet
%                         xs = xb(2) - sqrt( 2*hc*dx.^2/db );
%                         hmax = db / dx * ( xb(2) - xs );
%                         fphys{2}(:, subid  , 1) = 0;
%                         fphys{2}(:, subid+1, 1) = cell.project_vert2node( [0, hmax]' );
%                 end
%                 fphys{2}(:, subid  , 2) = 0;
%                 fphys{2}(:, subid+1, 2) = 0;
%             end
%         end% function
    end
    
end

