classdef StdPrismTri < handle
    
    properties ( Constant )
        %> reference element type
        type = enumStdCell.PrismTri;
        %> num of vertices
        Nv = 6
        %> volume of reference elements
        LAV = 4
        %> vertices coordinates
        vr = [-1,  1, -1, -1,  1, -1]'
        vs = [-1, -1,  1, -1, -1,  1]'
        vt = [-1, -1, -1,  1,  1,  1]'
        %> num of vertices on each face
        Nfv = [4, 4, 4, 3, 3];
        %> index of interpolation points on each face
        FToV = [1,2,5,4; 2,3,6,5; 3,1,4,6; 1,2,3,0; 4,5,6,0]';
        %> num of faces
        Nface = 5
        %> element types of each face
        faceType = [enumStdCell.Quad, ...
            enumStdCell.Quad, ...
            enumStdCell.Quad, ...
            enumStdCell.Tri, ...
            enumStdCell.Tri ]
    end
    
    properties
        %> maximum order of basis function
        N
        %> maximum polynomial order in vertical direction
        Nz
    end
    
    properties ( SetAccess = protected )
        %> num of interpolation points (IP)
        Np, Nph, Npz
        %> coordinates of interpolation points
        r, s, t
        %> horizontal and vertical node coordinate
        r1, s1, t1
        %> node index of facial points
        Fmask
        %> number of facial interpolation points
        Nfp
        %> total number of face points
        TNfp
    end
    properties ( SetAccess = protected)
        %> Vandermonde matrix
        V, Vh, Vint
        %> project matrx from interpolation points to quadrature points
        Vq
        %> mass matrix
        M
        %> inverse of mass matrix
        invM
        %> derivative matrix with
        %> \f$ [Dr]_{ij} = \left.\frac{\partial l_j}{\partial r}\right|_{r_i} \f$
        Dr, Ds, Dt
    end
    
    properties ( SetAccess = protected )
        %> number of gauss quadrature points
        Nq
        %> coordinate of quadrature points
        rq, sq, tq
        %> integral weights for each quadrature points
        wq
    end
    
    methods
        % construction function of triangular prism reference element
        function obj = StdPrismTri(Nh, Nz)
            obj.N = Nh;
            obj.Nz = Nz;
            assembleNodeCoor( obj, Nh, Nz );
            
            [ obj.Nq, obj.rq, obj.sq, obj.tq, obj.wq ] ...
                = obj.evaluateQuadCoor( Nh, Nz );
            
            assembleVandMatrix( obj );
            [ obj.Vq ] = obj.assembleQuadratureMatrix( );
            [ obj.M, obj.invM ] = obj.assembleMassMatrix( );
            [ obj.Dr, obj.Ds, obj.Dt ] ...
                = obj.evaluateNodalDerivativeFunc(obj.r, obj.s, obj.t);
            
            obj.Nfp = zeros(obj.Nface, 1);
            for i = 1:3
                obj.Nfp(i) = ( Nh + 1 ) * ( Nz +1 );
            end
            obj.Nfp( [ 4, 5 ] ) = obj.Nph;
            
            [ obj.TNfp ] = sum(obj.Nfp);
            [ obj.Fmask ] = obj.assembleFacialNodeIndex();
        end
        
        [ fun ] = evaluateOrthogonalFunc(obj, N1, N2, ind, r, s, t);
        [ node_val ] = projectVert2Node(obj, vert_val);
        
        %> @brief Evaluate all the nodal basis function values at points
        %> @param[in] obj The StdCell class
        %> @param[in] r The node coordinate
        %> @param[in] s The node coordinate
        %> @param[in] t The node coordinate
        %> @param[out] func The basis function values at points
        function [ func ] = evaluateNodalFunc(obj, r, s, t)
            func = zeros(numel(r), obj.Np);
            for n = 1:obj.Np
                fh = obj.evaluateHorizontalOrthogonalFunc( obj.N, n, r, s );
                fv = obj.evaluateVerticalOrthogonalFunc( n, t );
                func(:, n) = fh .* fv;
            end
            func = func/obj.V;
        end
        
        %> @brief Evaluate the derivative nodal function values at points
        function [ fDr, fDs, fDt ] = evaluateNodalDerivativeFunc( obj, r, s, t )
            Nr = numel( r );
            Vr = zeros(Nr, obj.Np);
            Vs = zeros(Nr, obj.Np);
            Vt = zeros(Nr, obj.Np);
            for n = 1:obj.Np
                [Vr(:, n), Vs(:, n), Vt(:, n)] = ...
                    obj.evaluateDerivativeOrthogonalFunc(obj.N, obj.Nz, n, r, s, t);
            end
            fDr = Vr/obj.V;
            fDs = Vs/obj.V;
            fDt = Vt/obj.V;
        end
    end% methods
    
    methods ( Access=protected )
        [ dr, ds, dt ] = evaluateDerivativeOrthogonalFunc( obj, Nh, Nv, ind, r, s, t );
        [ Nq, rq, sq, tq, wq ] = evaluateQuadCoor( obj, Nh, Nv );
        [ f ] = evaluateHorizontalOrthogonalFunc( obj, N1, td, r, s );
        [ f ] = evaluateVerticalOrthogonalFunc( obj, td, t );
        
        %> @brief Assemble the interpolation matrix of Gauss quadrature nodes
        %> The elements of the quadratuer interpolation matrix is
        %> \f$ [V_q]_{i,j} = l_j(\xi_i) \f$
        %> where \f$ \xi_i \f$ is the ith Gauss quadrature nodes.
        function [ Vq ] = assembleQuadratureMatrix( obj )
            Vq = obj.evaluateNodalFunc( obj.rq, obj.sq, obj.tq );
        end
        
        %> @brief Assemble the Vandermonde matrix
        %> @details The Vandermonde matrix interpolate the orthgonal basis functions
        %> to the nodal basis functions, with
        %> \f$ [V]_{i,j} = P_j(r_i) \f$
        %> where \f$ r_i \f$ is the ith interpolation nodes and \f$ P_j \f$
        %> is the jth orthgonal function.
        assembleVandMatrix( obj )
        
        %> @brief Assemble the mass matrix
        function [ M, invM ] = assembleMassMatrix( obj )
            invV = inv(obj.V);
            M = (invV')*invV;
            invM = obj.V * obj.V';
        end% func
        
        Fmask = assembleFacialNodeIndex(obj)        
    end
    
end

