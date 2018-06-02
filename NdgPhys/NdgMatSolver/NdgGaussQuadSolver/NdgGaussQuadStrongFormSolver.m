%> @brief class of solving strong form equations with the Gauss quadrature formula
%> 
%> 
classdef NdgGaussQuadStrongFormSolver < handle
    properties
        %> inverse mass matrix of each element
        invM
        %> values of derivative basis function at quadrature points
        %> \f$$ [\mathcal{D}_r]_{i,j} = \left. \frac{\partial \varphi_j}{\partial r} \right|_{\mathbf{r}_i}  \f$$
        Dr, Ds, Dt
        
        %> the element in the vector equals to
        %> \f$ \left. \frac{\partial r}{\partial x} \right|_{\mathbf{x}_i} \cdot w_i \cdot J_i \f$
        %> where \f$ w_i \f$ and \f$ J_i \f$ are the volume integral weights and Jacobian 
        %> at each quadrature points
        rxwJ, rywJ, rzwJ

        %> the element in the vector equals to
        %> \f$ \left. \frac{\partial s}{\partial x} \right|_{\mathbf{x}_i} \cdot w_i \cdot J_i \f$
        %> where \f$ w_i \f$ and \f$ J_i \f$ are the volume integral weights and Jacobian 
        %> at each quadrature points
        sxwJ, sywJ, szwJ
        
        %> the element in the vector equals to
        %> \f$ \left. \frac{\partial t}{\partial x} \right|_{\mathbf{x}_i} \cdot w_i \cdot J_i \f$
        %> where \f$ w_i \f$ and \f$ J_i \f$ are the volume integral weights and Jacobian 
        %> at each quadrature points
        txwJ, tywJ, tzwJ
        
        %> values of basis functions at each quadrature points, where
        %> \f$ [ \mathcal{V}_q ]_{i,j} = \varphi_{ \mathbf{r}_j } \f$
        Vq
    end
    properties
        %> values of the basis functions at each edge quadrature points
        %> \f$$ [\mathcal{LIFT}_t]_{i,ej} = \varphi_i(\mathbf{r}_{ej}), \f$$
        %> where \f$ \mathbf{r}_{ej} \f$ is the location of quadrature points on edges.
        LIFT
        %> outward normal vector at the surface integral quadrature points
        nx, ny, nz
        
        %> the element in the vector equals to
        %> \f$ w_{ei} \cdot Js_{ei} \f$
        %> where \f$ w_{ei} \f$ and \f$ Js_{ei} \f$ are the surface integral weights and Jacobian 
        %> at each quadrature points on edges
        wJs
        %> project nodal values to face quadrature points
        Vfq
        %> project face values to face quadrature points
        FVfq
        %> total number of quadrature points for surface integration
        TNfq
    end
    
    methods
        function obj = NdgGaussQuadStrongFormSolver( phys )
            
            obj.TNfq = cell( phys.Nmesh, 1 );
            for m = 1:phys.Nmesh
                mesh = phys.meshUnion( m );
                stdcell = mesh.cell;
                % count the total quadrature points on each edge
                TNfq = 0;
                for f = 1:stdcell.Nface
                    fcell = getStdCell( stdcell.N, stdcell.faceType(f) );
                    TNfq = TNfq + fcell.Nq;
                end
                [ obj.TNfq{m} ] = TNfq;
                [ obj.Vq{m} ] = stdcell.Vq';
                [ obj.invM{m} ] = obj.assembleInverseMassMatrix( mesh );
                [ obj.Vfq{m} ] = obj.assembleVandMatrixFaceQuadrature( mesh.cell, TNfq );
                [ obj.FVfq{m} ] = obj.assembleFacialVandMatrixFaceQuadrature( mesh.cell, TNfq );
                [ obj.Dr{m}, obj.Ds{m}, obj.Dt{m} ] = obj.assembleDerivativeMatrix( mesh );                
                [ nx, ny, nz, Js ] = obj.assembleNormalVector( mesh );
                % project the outward normal vector from edge nodes to quadrature nodes.
                obj.nx{m} = obj.FVfq{m} * nx;
                obj.ny{m} = obj.FVfq{m} * ny;
                obj.nz{m} = obj.FVfq{m} * nz;
                [ obj.wJs{m} ] = obj.assembleFaceQuadratureWeight( mesh, TNfq, obj.FVfq{m}*Js );
                [ obj.rxwJ{m}, obj.rywJ{m}, obj.rzwJ{m}, ...
                    obj.sxwJ{m}, obj.sywJ{m}, obj.szwJ{m}, ...
                    obj.txwJ{m}, obj.tywJ{m}, obj.tzwJ{m} ] = obj.assembleJacobianFactor( mesh );
                
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
        
        function [ rxwJ, rywJ, rzwJ, sxwJ, sywJ, szwJ, txwJ, tywJ, tzwJ ] ...
                = assembleJacobianFactor( mesh )
            
            [ rx, ry, rz, sx, sy, sz, tx, ty, tz, J ] ...
                = mesh.cell.assembleJacobianMatrix( mesh.x, mesh.y, mesh.z );
            [ Jq ] = mesh.cell.project_node2quad( J );
            [ wJ ] = bsxfun(@times, mesh.cell.wq, Jq);
            
            rxwJ = wJ.*( mesh.cell.project_node2quad( rx ) );
            rywJ = wJ.*( mesh.cell.project_node2quad( ry ) );
            rzwJ = wJ.*( mesh.cell.project_node2quad( rz ) );
            txwJ = wJ.*( mesh.cell.project_node2quad( tx ) );
            tywJ = wJ.*( mesh.cell.project_node2quad( ty ) );
            tzwJ = wJ.*( mesh.cell.project_node2quad( tz ) );
            sxwJ = wJ.*( mesh.cell.project_node2quad( sx ) );
            sywJ = wJ.*( mesh.cell.project_node2quad( sy ) );
            szwJ = wJ.*( mesh.cell.project_node2quad( sz ) );
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

