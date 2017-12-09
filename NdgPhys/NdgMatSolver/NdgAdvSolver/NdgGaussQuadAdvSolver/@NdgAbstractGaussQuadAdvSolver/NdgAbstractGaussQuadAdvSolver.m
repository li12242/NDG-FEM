classdef NdgAbstractGaussQuadAdvSolver < NdgAbstractAdvSolver
    
    properties
        invM
        Dr
        Ds
        Dt
        LIFT
        nx
        ny
        nz
        rxwJ
        rywJ
        rzwJ
        sxwJ
        sywJ
        szwJ
        txwJ
        tywJ
        tzwJ
        JswJ
        wJs
        %> project nodal values to face quadrature points
        Vfq
        %> project nodal values to quadrature points
        Vq
        %> project face values to face quadrature points
        FVfq
    end
    
    methods
        function obj = NdgAbstractGaussQuadAdvSolver( phys )
            obj = obj@NdgAbstractAdvSolver( phys );
            
            for m = 1:phys.Nmesh
                mesh = phys.meshUnion( m );
                cell = mesh.cell;
                % count the total quadrature points on each face
                TNfq = 0;
                for f = 1:cell.Nface
                    fcell = getStdCell( cell.N, cell.faceType(f) );
                    TNfq = TNfq + fcell.Nq;
                end
                
                [ obj.Vq{m} ] = cell.Vq;
                [ obj.invM{m} ] = obj.assembleInverseMassMatrix( mesh );
                [ obj.Vfq{m} ] = obj.assembleVandMatrixFaceQuadrature( mesh.cell, TNfq );
                [ obj.FVfq{m} ] = obj.assembleFacialVandMatrixFaceQuadrature( mesh.cell, TNfq );
                [ obj.Dr{m}, obj.Ds{m}, obj.Dt{m} ] = obj.assembleDerivativeMatrix( mesh );
                
                [ rx, ry, rz, sx, sy, sz, tx, ty, tz, J ] = obj.assembleJacobianFactor( mesh );
                [ nx, ny, nz, Js ] = obj.assembleNormalVector( mesh );
                
                obj.nx{m} = obj.FVfq{m} * nx;
                obj.ny{m} = obj.FVfq{m} * ny;
                obj.nz{m} = obj.FVfq{m} * nz;
                [ obj.wJs{m} ] = obj.assembleFaceQuadratureWeight( mesh, TNfq, Js );
                [ Jq ] = cell.project_node2quad( J );
                [ wJ ] = bsxfun(@times, mesh.cell.wq, Jq);
                obj.rxwJ{m} = wJ.*( mesh.cell.project_node2quad( rx ) );
                obj.rywJ{m} = wJ.*( mesh.cell.project_node2quad( ry ) );
                obj.rzwJ{m} = wJ.*( mesh.cell.project_node2quad( rz ) );
                obj.txwJ{m} = wJ.*( mesh.cell.project_node2quad( tx ) );
                obj.tywJ{m} = wJ.*( mesh.cell.project_node2quad( ty ) );
                obj.tzwJ{m} = wJ.*( mesh.cell.project_node2quad( tz ) );
                obj.sxwJ{m} = wJ.*( mesh.cell.project_node2quad( sx ) );
                obj.sywJ{m} = wJ.*( mesh.cell.project_node2quad( sy ) );
                obj.szwJ{m} = wJ.*( mesh.cell.project_node2quad( sz ) );
                
                [ obj.LIFT{m} ] = obj.assembleLiftMatrix( mesh, TNfq );
            end
        end
    end
    
    methods( Static )
        function [ wJs ] = assembleFaceQuadratureWeight( mesh, TNfq, Js )
            cell = mesh.cell;
            sk = 1;
            ws = zeros( TNfq, 1 );
            for f = 1:cell.Nface
                fcell = getStdCell( cell.N, cell.faceType(f) );
                ws(sk:(sk+fcell.Nq-1) ) = fcell.wq;
                sk = sk + fcell.Nq;
            end
            wJs = bsxfun(@times, ws, Js);
        end
        
        function [ rx, ry, rz, sx, sy, sz, tx, ty, tz, J ] = assembleJacobianFactor( mesh )
            [ rx, ry, rz, sx, sy, sz, tx, ty, tz, J ] ...
                = mesh.cell.assembleJacobianMatrix( mesh.x, mesh.y, mesh.z );
        end
        
        function [ nx, ny, nz, Js ] = assembleNormalVector( mesh )
            [ nx, ny, nz, Js ] = mesh.cell.assembleNormalVector( mesh.x, mesh.y, mesh.z );
        end
        
        function [ Vfq ] = assembleVandMatrixFaceQuadrature( cell, TNfq )
            sk = 1;
            rfq = zeros( TNfq, 1 );
            sfq = zeros( TNfq, 1 );
            tfq = zeros( TNfq, 1 );
            for f = 1:cell.Nface
                fcell = getStdCell( cell.N, cell.faceType(f) );
                vr = cell.vr( cell.FToV(:, f) );
                vs = cell.vs( cell.FToV(:, f) );
                vt = cell.vt( cell.FToV(:, f) );
                rq = fcell.project_vert2quad( vr );
                sq = fcell.project_vert2quad( vs );
                tq = fcell.project_vert2quad( vt );
                rfq( sk:(sk+fcell.Nq-1) ) = rq(:);
                sfq( sk:(sk+fcell.Nq-1) ) = sq(:);
                tfq( sk:(sk+fcell.Nq-1) ) = tq(:);
                sk = sk+fcell.Nq;
            end
            
            Vfq = cell.nodal_func( rfq, sfq, tfq );
        end
        
        function [ FVfq ] = assembleFacialVandMatrixFaceQuadrature( cell, TNfq )
            sk = 1;
            sp = 1;
            FVfq = zeros( TNfq, cell.TNfp );
            for f = 1:cell.Nface
                fcell = getStdCell( cell.N, cell.faceType(f) );
                FVfq( sk:(sk+fcell.Nq-1), sp:(sp+fcell.Np-1) ) = fcell.Vq;
                sk = sk+fcell.Nq;
                sp = sp+fcell.Np;
            end
        end
        
        function [ invM ] = assembleInverseMassMatrix( mesh )
            cell = mesh.cell;
            Np = cell.Np;
            K = mesh.K;
            invM = zeros( Np, Np, K );
            for k = 1:K
                Jq = cell.project_node2quad( mesh.J(:, k) );
                M = ( cell.Vq' * diag( Jq.*cell.wq ) ) * cell.Vq;
                invM(:, :, k) = inv( M );
            end
        end
        
        function [ Dr, Ds, Dt ] = assembleDerivativeMatrix( mesh )
            cell = mesh.cell;
            [ Dr, Ds, Dt ] = cell.nodal_derivative_func( cell.rq, cell.sq, cell.tq );
        end
        
        function [ LIFT ] = assembleLiftMatrix( mesh, TNfq )
            cell = mesh.cell;
            sk = 1;
            LIFT = zeros( cell.Np, TNfq );
            for f = 1:cell.Nface
                fcell = getStdCell( cell.N, cell.faceType(f) );
                vr = cell.vr( cell.FToV(:, f) );
                vs = cell.vs( cell.FToV(:, f) );
                vt = cell.vt( cell.FToV(:, f) );
                rq = fcell.project_vert2quad( vr );
                sq = fcell.project_vert2quad( vs );
                tq = fcell.project_vert2quad( vt );
                LIFT(:, sk:(sk+fcell.Nq-1) ) = ( cell.nodal_func( rq, sq, tq ) )';
                sk = sk + fcell.Nq;
            end
        end
    end
    
end

