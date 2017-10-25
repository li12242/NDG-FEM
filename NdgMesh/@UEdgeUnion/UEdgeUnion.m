%> @brief Edge information on unstructured meshes.
%
%> And here we can put some more detailed informations about the class.
% ======================================================================
%> This class is part of the NDG-FEM software.
%> @author li12242, Tianjin University, li12242@tju.edu.cn
% ======================================================================
classdef UEdgeUnion < handle
    
    properties(SetAccess=protected)
        %> pointer to unstructured mesh object.
        umesh
        %> pointer to edge cell object
        bcell
        %> number of edges
        M
        %> vertex index on each edge
        FToV
        %> local and adjacent cell index
        FToE
        %> local face index of local and adjacent cell
        FToF
        %> mesh index of local and adjacent cell
        FToM
        %> boundary condition types of edges
        edgeType
        
        %> project the boundary nodal values to the quadrature points for local cell
        NToQ1
        %> project the boundary nodal values to the quadrature points for adjacent cell
        NToQ2
        
        %> determination of Jacobian matrixs
        Jq
        %> x component of the outward normal vector
        nxq
        %> y component of the outward normal vector
        nyq
        %> z component of the outward normal vector
        nzq
    end
    
    methods(Abstract, Hidden=true, Access=protected)
        [ bcell ] = setStdEdgeCell( obj, N );
        [ nxq, nyq, nzq, Jq ] = assembleQuadPointScale( obj );
    end
    
    methods
        function obj = UEdgeUnion( mesh, meshId, BCToV)
            obj = obj.constructor( mesh, meshId, BCToV );
        end
    end
    
    methods(Hidden, Access=protected)
        function obj = constructor(obj, mesh, meshId, BCToV)
            [ obj.umesh ] = mesh;
            [ obj.bcell ] = setStdEdgeCell( obj );
            [ obj.M, obj.FToM, obj.FToE, obj.FToF, obj.FToV, obj.edgeType ] ...
                = assembleEdgeConnect( obj, meshId, BCToV );
            [ obj.NToQ1, obj.NToQ2 ] = assembleNodeProject( obj );
            [ obj.nxq, obj.nyq, obj.nzq, obj.Jq ] ...
                = assembleQuadPointScale( obj );
        end
        
        [ Nedge, FToM, FToE, FToF, FToV, edgeType ] ...
            = assembleEdgeConnect( obj, meshId, BCToV );
        
        function [ NToQ1, NToQ2 ] = assembleNodeProject( obj )
            mesh = obj.umesh;
            m1 = obj.FToM(1);
            m2 = obj.FToM(2);
            cell1 = mesh( m1 ).cell;
            cell2 = mesh( m2 ).cell;
            Nq = obj.bcell.Nq;
            NToQ1 = zeros( Nq, cell1.Np, obj.M );
            NToQ2 = zeros( Nq, cell2.Np, obj.M );

            for n = 1:obj.M
                k1 = obj.FToE(1, n);
                k2 = obj.FToE(2, n);
                v1 = mesh( m1 ).EToV(:, k1);
                v2 = mesh( m2 ).EToV(:, k2);
                f1 = obj.FToF(1, n);
                f2 = obj.FToF(2, n);
                vert1 = v1( cell1.FToV(:, f1) );
                try
                vert2 = v2( cell2.FToV(:, f2) );
                catch
                    keyboard
                end
                
                [ vert, vertId1 ] = sort(vert1);
                [ ~, vertId2 ] = sort(vert2);
                obj.FToV(:, n) = vert;
                                
                rv1 = cell1.vr( cell1.FToV(vertId1, f1) );
                rv2 = cell2.vr( cell2.FToV(vertId2, f2) );
                sv1 = cell1.vs( cell1.FToV(vertId1, f1) );
                sv2 = cell2.vs( cell2.FToV(vertId2, f2) );
                rq1 = obj.bcell.project_vert2quad(rv1);
                rq2 = obj.bcell.project_vert2quad(rv2);
                sq1 = obj.bcell.project_vert2quad(sv1);
                sq2 = obj.bcell.project_vert2quad(sv2);
                for i = 1:cell1.Np
                    NToQ1(:,i,n) = cell1.orthogonal_func(cell1.N, i, rq1, sq1, 0);
                end
                for i = 1:cell2.Np
                    NToQ2(:,i,n) = cell2.orthogonal_func(cell2.N, i, rq2, sq2, 0);
                end
                NToQ1(:,:,n) = NToQ1(:,:,n)/cell1.V;
                NToQ2(:,:,n) = NToQ2(:,:,n)/cell2.V;
            end% for
        end
    end
end

