%> @brief standard cell class.
%
%> And here we can put some more detailed informations about the class.
% ======================================================================
%> This class is part of the NDGOM software. 
%> @author li12242, Tianjin University, li12242@tju.edu.cn
% ======================================================================
classdef StdCell < handle
    properties
        %> order of basis function
        N
    end
    
    properties(Constant, Abstract)
        %> standard cell type
        type
        %> number of vertex
        Nv
        %> vertex coordinate
        vr
        %> vertex coordinate
        vs
        %> vertex coordinate
        vt
        %> length/area/volume of standard cell
        vol
        %> number of vertex on each face
        Nfv
        %> veretx list on each face
        FToV
        %> number of faces
        Nface
        %> standard cell types of each face
        faceType
    end
    
    properties(SetAccess = protected)
        %> number of interpolation points (IP)
        Np
        %> coordinates of interpolation points
        r
        %> coordinates of interpolation points
        s
        %> coordinates of interpolation points
        t
        %> node index of facial points
        Fmask
        %> number of facial interpolation points
        Nfp
        %> total number of face points
        TNfp
    end
    properties(Hidden, SetAccess = protected)
        %> Vandermonde matrix
        V
        %> mass matrix
        M
        %> derivative matrix, 
        %> \f$ [Dr]_{ij} = \left.\frac{\partial l_j}{\partial r}\right|_{r_i} \f$
        Dr
        %> derivative matrix, 
        %> \f$ [Ds]_{ij} = \left.\frac{\partial l_j}{\partial s}\right|_{r_i} \f$
        Ds
        %> derivative matrix, 
        %> \f$ [Dt]_{ij} = \left.\frac{\partial l_j}{\partial t}\right|_{r_i} \f$
        Dt
        %> lift matrix, \f$ LIFT = M^{-1} \cdot M_e \f$
        LIFT
        %> project matrx from interpolation points to quadrature points
        Vq
    end

    properties(Hidden, SetAccess = protected)
        %> number of gauss quadrature points
        Nq
        %> coordinate of quadrature points
        rq
        %> coordinate of quadrature points
        sq
        %> coordinate of quadrature points
        tq
        %> integral weights for each quadrature points
        wq
    end
    
    methods(Abstract, Access=protected)
        %> get the total number and coordinate of interpolation points
        [ Np,r,s,t ] = node_coor_func(obj, N)
        %> get the total number and coordinate of gauss quadrature points
        [ Nq,rq,sq,tq ] = quad_coor_func(obj, N)
        %> get the derivative of basis function at each interpolation points
        [dr, ds, dt] = derivative_orthogonal_func(obj, N, ind, r, s, t)
    end
    
    methods(Abstract)
        fun = orthogonal_func(obj, N, ind, r, s, t);
        node_val = project_vert2node(obj, vert_val);
    end
    
    methods
        function obj = StdCell(N)
            obj.N = N;
            [ obj.Np, obj.r, obj.s, obj.t ] = obj.node_coor_func(N);
            [ obj.V ] = obj.Vandmode_Matrix( @obj.orthogonal_func );
            [ obj.M ] = obj.Mass_Matrix();
            [ obj.Dr, obj.Ds, obj.Dt ] = ...
                obj.Derivative_Matrix( @obj.derivative_orthogonal_func );
            
            [ obj.Nq, obj.rq, obj.sq, obj.tq, obj.wq ] = quad_coor_func(obj, N);
            [ obj.Vq ] = quad_matrix( obj );
            % face
            if obj.Nface > 0
                obj.Nfp = zeros(obj.Nface, 1);
            else
                obj.Nfp = 1;
            end
            for f = 1:obj.Nface
                cell = getStdCell(obj.N, obj.faceType(f));
                obj.Nfp(f) = cell.Np;
            end
            obj.TNfp = sum(obj.Nfp);
            obj.Fmask = obj.Face_node();
            obj.LIFT = obj.Lift_Matrix();
        end

        %> project the node values on the quadrature nodes
        function quad_val = project_node2quad(obj, node_val)
            quad_val = obj.Vq * node_val;
        end% func

        function quad_val = project_vert2quad(obj, vert_val)
            node_val = obj.project_vert2node(vert_val);
            quad_val = obj.project_node2quad(node_val);
        end
    end% methods
   
    methods(Hidden, Access = protected)
        function [ Vq ] = quad_matrix( obj )
            Vq = zeros(obj.Nq, obj.Np);
            for n = 1:obj.Np
                Vq(:, n) = obj.orthogonal_func(obj.N, n, obj.rq, obj.sq, obj.tq);
            end
            Vq = Vq/obj.V;
        end
        
        function V = Vandmode_Matrix(obj, orthogonal_func)
            V = zeros(obj.Np, obj.Np);
            for n = 1:obj.Np
                V(:, n) = orthogonal_func(obj.N, n, obj.r, obj.s, obj.t);
            end% for
        end% func
        
        function M = Mass_Matrix( obj )
            invV = inv(obj.V);
            M = (invV')*invV;
        end% func
        
        function [Dr, Ds, Dt] = Derivative_Matrix(obj, deri_orthogonal_func)
            Vr = zeros(obj.Np, obj.Np);
            Vs = zeros(obj.Np, obj.Np);
            Vt = zeros(obj.Np, obj.Np);
            for n = 1:obj.Np
                [Vr(:, n), Vs(:, n), Vt(:, n)] = deri_orthogonal_func...
                    (obj.N, n, obj.r, obj.s, obj.t);
            end
            Dr = Vr/obj.V; 
            Ds = Vs/obj.V; 
            Dt = Vt/obj.V;
        end% func
        
        function Fmask = Face_node(obj)
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
        
        function LIFT = Lift_Matrix(obj)
            Mes = zeros(obj.Np, obj.TNfp);
            sk = 1;
            for f = 1:obj.Nface
                cell = getStdCell(obj.N, obj.faceType(f));
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