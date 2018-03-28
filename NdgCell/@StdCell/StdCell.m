%> @brief Abstract class of the standard cells.
%>
%> The basic properties of the StdCell object include its geometry 
%> information and the matrix coefficients.
%>
%> All StdCell use the nodal interpolation polynomial as the basis
%> functions, while the LGL nodes are choosen as the interpolation nodes.
%> To create a new StdCell object, the user has to input the maximum degree
%> of the basis functions, such as
%> @code
%>  tri = StdTri(2); % the maximum degree of the basis polynomial is 2.
%> @endcode
%> 
% ======================================================================
%> The public interface of StdCell includes:
%> @code
%>   [ func ] = orthogonal_func(obj, N, ind, r, s, t); //evaluate the orthogonal function values at points
%>   [ func ] = nodal_func(obj, r, s, t); // evaluate all the nodal basis function values at points
%>   [ Dx, Dy, Dz ] = nodal_deri_func(obj, x, y, z);
%>   [ node_val ] = project_vert2node(obj, vert_val); // evaluate the node values from the vertice values
%>   [ quad_val ] = project_node2quad(obj, node_val); // evaluate the quadrature node values from the node values
%>   [ quad_val ] = project_vert2quad(obj, vert_val); // evaluate the quadrature node values from the vertice values
%> @endcode
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
        LAV
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
    properties( SetAccess = protected)
        %> Vandermonde matrix
        V
        %> project matrx from interpolation points to quadrature points
        Vq
        %> mass matrix
        M
        %> inverse of mass matrix
        invM
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
    end

    properties( SetAccess = protected )
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
        [ Nq,rq,sq,tq,wq ] = quad_coor_func(obj, N)
        %> get the derivative of orthogonal function at each interpolation points
        [dr, ds, dt] = derivative_orthogonal_func(obj, N, ind, r, s, t);
    end
    
    methods(Abstract)
        %> @brief Return the value of the orthogonal function at nodes (r, s, t)
        %> @param[in] obj StdCell class
        %> @param[in] N The maximum number of the basis function
        %> @param[in] ind The index of the orthgonal function
        %> @param[in] r Coordinate of the nodes
        %> @param[in] s Coordinate of the nodes
        %> @param[in] t Coordinate of the nodes
        [ fun ] = orthogonal_func(obj, N, ind, r, s, t);
        %> @brief Project the scalar field from the cell vertices to the interpolation nodes.
        %> @param[in] obj StdCell class
        %> @param[in] vert_val the vertice values, the first dimension of
        %> this variable should be equal to the Nv of the object.
        [ node_val ] = project_vert2node(obj, vert_val);
        %> @brief Calculate the 
        %> @param[in]
        assembleJacobianMatrix( obj, x, y, z );
        %> @brief Assemble the outword normal vectors.
        %> @param[in]
        %> @param[in]
        assembleNormalVector( obj, x, y, z )
    end
    
    methods
        %> Construction function of the StdCell class
        %> @param[in] N The maximum degree of the basis functions
        function obj = StdCell(N)
            [ obj.N ] = N;
            [ obj.Np, obj.r, obj.s, obj.t ] = obj.node_coor_func( N );            
            [ obj.Nq, obj.rq, obj.sq, obj.tq, obj.wq ] = obj.quad_coor_func( N );
            [ obj.V ] = obj.assembleVandMatrix( @obj.orthogonal_func );
            [ obj.Vq ] = obj.assembleQuadratureMatrix();
            [ obj.M, obj.invM ] = obj.assembleMassMatrix();
            [ obj.Dr, obj.Ds, obj.Dt ] = obj.nodal_derivative_func(obj.r, obj.s, obj.t);
            %[ obj.Drq, obj.Dsq, obj.Dtq ] ...
            %    = obj.assembleQuadratureDerivativeMatrix( @obj.derivative_orthogonal_func );
            
            % get the number of nodes on each face
            if obj.Nface > 0
                obj.Nfp = zeros(obj.Nface, 1);
            else
                obj.Nfp = 1;
            end
            for f = 1:obj.Nface
                cell = getStdCell(obj.N, obj.faceType(f));
                obj.Nfp(f) = cell.Np;
            end
            [ obj.TNfp ] = sum(obj.Nfp);
            [ obj.Fmask ] = obj.assembleFacialNodeIndex();
            [ obj.LIFT ] = obj.assembleLiftMatrix();
        end
        
        %> @brief Evaluate all the nodal basis function values at points
        %> @param[in] obj The StdCell class
        %> @param[in] r The node coordinate
        %> @param[in] s The node coordinate
        %> @param[in] t The node coordinate
        %> @param[out] func The basis function values at points
        function [ func ] = nodal_func(obj, r, s, t)
            func = zeros(numel(r), obj.Np);
            for n = 1:obj.Np
                func(:, n) = obj.orthogonal_func(obj.N, n, r, s, t);
            end
            func = func/obj.V;
        end
        
        function [ dfr, dfs, dft ] = orthogonal_derivative_func(obj, ind, r, s, t)
            [ dfr, dfs, dft ] = obj.derivative_orthogonal_func( obj.N, ind, r, s, t );
        end
        
        %> @brief Evaluate the derivative nodal function values at points
        function [ fDr, fDs, fDt ] = nodal_derivative_func( obj, r, s, t )
            Nr = numel( r );
            Vr = zeros(Nr, obj.Np);
            Vs = zeros(Nr, obj.Np);
            Vt = zeros(Nr, obj.Np);
            for n = 1:obj.Np
                [Vr(:, n), Vs(:, n), Vt(:, n)] = ...
                    obj.derivative_orthogonal_func(obj.N, n, r, s, t);
            end
            fDr = Vr/obj.V; 
            fDs = Vs/obj.V; 
            fDt = Vt/obj.V;
        end

%> @brief Project the scalar field from the interpolation nodes to the Gauss quadrature nodes
%> @param[in] obj The StdCell class
%> @param[in] node_val The values on these interpolation nodes
        function quad_val = project_node2quad(obj, node_val)
            quad_val = obj.Vq * node_val;
        end% func

%> @brief Project the scalar field from the vertices to the Gauss quadrature nodes
%> @param[in] obj The StdCell class
%> @param[in] vert_val The values on these vertices
        function quad_val = project_vert2quad(obj, vert_val)
            node_val = obj.project_vert2node(vert_val);
            quad_val = obj.project_node2quad(node_val);
        end
    end% methods
   
    methods(Hidden, Access = protected)
        %> @brief Assemble the interpolation matrix of Gauss quadrature nodes
        %>
        %> The elements of the quadratuer interpolation matrix is 
        %> \f$ [V_q]_{i,j} = l_j(\xi_i) \f$
        %> where \f$ \xi_i \f$ is the ith Gauss quadrature nodes.
        function [ Vq ] = assembleQuadratureMatrix( obj )
            Vq = obj.nodal_func( obj.rq, obj.sq, obj.tq );
        end
        
        %> @brief Assemble the Vandermonde matrix
        %> 
        %> The Vandermonde matrix interpolate the orthgonal basis functions
        %> to the nodal basis functions, with 
        %> \f$ [V]_{i,j} = P_j(r_i) \f$
        %> where \f$ r_i \f$ is the ith interpolation nodes and \f$ P_j \f$
        %> is the jth orthgonal function.
        function V = assembleVandMatrix(obj, orthogonal_func)
            V = zeros(obj.Np, obj.Np);
            for n = 1:obj.Np
                V(:, n) = orthogonal_func(obj.N, n, obj.r, obj.s, obj.t);
            end% for
        end% func
        
        %> @brief Assemble the mass matrix
        %> 
        function [ M, invM ] = assembleMassMatrix( obj )
%             [ Nq, rq, sq, tq, wq ] = quad_coor_func( obj, 2*obj.N );
%             [ fq ] = nodal_func( obj, rq, sq, tq );
%             M = zeros( obj.Np, obj.Np );
%             for n = 1:obj.Np
%                 for m = 1:obj.Np
%                     M(n, m) = sum( wq(:) .* fq(:, m) .* fq(:, n) );
%                 end
%             end
%             invM = inv(M);
            invV = inv(obj.V);
            M = (invV')*invV;
            invM = obj.V * obj.V';
        end% func
        
        function Fmask = assembleFacialNodeIndex(obj)
            maxnfp = max(obj.Nfp);
            Fmask = zeros(maxnfp, obj.Nface);
            for f = 1:obj.Nface
                nfv = obj.Nfv(f);
                % get vertex index on face f
                rv = obj.vr( obj.FToV(1:nfv, f) );
                sv = obj.vs( obj.FToV(1:nfv, f) );
                tv = obj.vt( obj.FToV(1:nfv, f) );
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
                fcell = getStdCell(obj.N, obj.faceType(f));
                row = obj.Fmask(:, f);
                row = row(row ~= 0);
                for n = 1:fcell.Np
                    Mes(row, sk) = fcell.M(:, n);
                    sk = sk + 1;
                end
            end
            LIFT = obj.invM * Mes;
        end
    end% methods
end