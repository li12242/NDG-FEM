%> @brief Super class for standard cell classes.
%
%> And here we can put some more detailed informations about the class.
% ======================================================================
%> This class is part of the NDG-FEM software. 
%> @author li12242, Tianjin University
%> @email li12242@tju.edu.cn
% ======================================================================
classdef StdCell < handle
    
    properties
        %> The order of the cell
        N
    end
    
    properties(Constant, Abstract)
        %> cell type
        type
        %> number of vertices
        Nv
        %> vertices coordinate
        vr, vs, vt
        %> length/area/volume of the standard cell
        vol
        %> number of vertices on each face
        Nfv
        %> vertices index on each face
        FToV
        %> number of faces
        Nface
        %> cell types of each face
        faceType
    end
    
    properties(SetAccess = protected)
        %> number of interpolation points
        Np
        %> coordinates of interpolation points
        r, s, t
        %> vandermonde matrix
        V
        %> mass matrix
        M
        %> derivative matrix
        Dr, Ds, Dt        
    end
    
    properties(SetAccess = protected)
        %> number of quadrature
        Nq
        %> quadrature points
        rq, sq, tq
        %> the basis function vales at the quadrature points
        Vq
        %> integral weights for each IPS
        wq
    end
    
    properties(SetAccess = protected)
        %> index of the IPS on the edges
        Fmask
        %> number of IPS on each edge
        Nfp
        %> total number of the edge IPPS
        TNfp
        %> LIFT matrix
        LIFT
    end
    
    methods(Abstract, Access=protected)
        %> set the number of interpolation nodes and their coordinates
        [Np, r,s,t] = node_coor_func(obj, N)
        %> calculate the values of the derivative orthogonal function
        [dr, ds, dt] = derivative_orthogonal_func(obj, N, ind, r, s, t)
        %> collect quadrature nodes coordinate and number
        [Nq, rq, sq, tq, wq] = quadrature_node_func(obj, N)
    end
    
    methods(Abstract, Access=public)
        %> calculate the value of the orthogonal function
        fun = orthogonal_func(obj, N, ind, r, s, t);
        %> project the vertex values to the node values
        node_val = project_vert2node(obj, vert_val);        
    end
    
    methods
        %> Construct the object with specific order N.
        function obj = StdCell(N)
            obj.N = N;
            % volume
            [ obj.Np, obj.r, obj.s, obj.t ] = node_coor_func(obj, N);
            [ obj.V ] = assembleVandmodeMatrix(obj, @obj.orthogonal_func);
            [ obj.M ] = assembleMassMatrix( obj );
            [ obj.Dr, obj.Ds, obj.Dt ] = ...
                assembleDerivativeMatrix(obj, @obj.derivative_orthogonal_func);
            
            % face
            if obj.Nface > 0
                obj.Nfp = zeros(obj.Nface, 1);
                cell = getStdCell( repmat(obj.N, obj.Nface, 1), obj.faceType );
                obj.Nfp = [ cell(:).Np ];
            else
                obj.Nfp = 1;
            end
            
            [ obj.Nq, obj.rq, obj.sq, obj.tq, obj.wq ] = quadrature_node_func(obj, N);
            [ obj.Vq ] = assembleQuadratureMapMatrix(obj, @obj.orthogonal_func);
            
            obj.TNfp = sum( obj.Nfp );
            obj.Fmask = assembleFaceMask( obj );
            obj.LIFT = assembleLiftMatrix( obj );
        end
        
        %> project the node values on the quadrature nodes
        function quad_val = project_node2quad(obj, node_val)
            quad_val = obj.Vq * node_val;
        end% func
    end
    
    methods(Hidden, Access = protected)
        %> assemble the map matrix from the IPPS to the quadrature nodes
        function Vq = assembleQuadratureMapMatrix(obj, orthogonal_func)
            Vq = zeros(obj.Nq, obj.Np);
            for n = 1:obj.Np
                Vq(:, n) = orthogonal_func(obj.N, n, obj.rq, obj.sq, obj.tq);
            end
            Vq = Vq/obj.V;
        end
        
        function V = assembleVandmodeMatrix(obj, orthogonal_func)
            V = zeros(obj.Np, obj.Np);
            for n = 1:obj.Np
                V(:, n) = orthogonal_func(obj.N, n, obj.r, obj.s, obj.t);
            end% for
        end% func
        
        function M = assembleMassMatrix(obj)
            invV = inv(obj.V);
            M = (invV')*invV;
        end
        
        %> assemble the derivative matrix
        function [Dr, Ds, Dt] = assembleDerivativeMatrix(...
                obj, deri_orthogonal_func)
            
            Vr = zeros(obj.Np, obj.Np);
            Vs = zeros(obj.Np, obj.Np);
            Vt = zeros(obj.Np, obj.Np);
            for n = 1:obj.Np
                [Vr(:, n), Vs(:, n), Vt(:, n)] = deri_orthogonal_func...
                    (obj.N, n, obj.r, obj.s, obj.t);
            end
            Dr = Vr/obj.V; Ds = Vs/obj.V; Dt = Vt/obj.V;
        end% func
        
        %> assemble the Fmask matrix
        function Fmask = assembleFaceMask(obj)
            maxnfp = max(obj.Nfp);
            Fmask = zeros(maxnfp, obj.Nface);
            for f = 1:obj.Nface
                nfv = obj.Nfv(f);
                % get vertex index on face f
                rv = obj.vr(obj.FToV(1:nfv, f));
                sv = obj.vs(obj.FToV(1:nfv, f));
                tv = obj.vt(obj.FToV(1:nfv, f));
                if(isrow(rv)) rv = rv'; end
                if(isrow(sv)) sv = sv'; end
                if(isrow(tv)) tv = tv'; end
                % get the nodes on face f
                cell = getStdCell(obj.N, obj.faceType(f));
                fr = cell.project_vert2node(rv);
                fs = cell.project_vert2node(sv);
                ft = cell.project_vert2node(tv);
                % get the nodes index
                for n = 1:obj.Nfp(f)
                    dis = (fr(n) - obj.r).^2 + (fs(n) - obj.s).^2 + (ft(n) - obj.t).^2;
                    ind = find(dis < 1e-10);
                    Fmask(n, f) = ind;
                end
            end
        end% func
        
        function LIFT = assembleLiftMatrix(obj)
            Mes = zeros(obj.Np, obj.TNfp);
            sk = 1;
            for f = 1:obj.Nface
                cell = ndg_lib.get_std_cell(obj.N, obj.faceType(f));
                row = obj.Fmask(:, f);
                row = row(row ~= 0);
                for n = 1:cell.Np
                    try
                        Mes(row, sk) = cell.M(:, n);
                    catch
                        keyboard
                    end
                    sk = sk + 1;
                end
            end
            LIFT = (obj.V*(obj.V)')*Mes;
        end
    end
end

