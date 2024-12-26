classdef refelement
    % refelement reference element for the nodal discontinuous galerkin method
    %
    % :param N: polynomial order

    properties
        % maximum order of basis function
        numOrder
    end

    properties ( Constant, Abstract )
        % num of vertice
        nVert
        % number of faces
        nFace
        % reference cell vertice coordinate - r
        refVR
        % reference cell vertice coordinate - s
        refVS
        % reference cell vertice coordinate - t
        refVT
        % length/area/volume of reference cell
        refSize
        % num of vertice on each face
        nFaceVert
        % veretx list on each face
        FToV
        % standard cell types of each face
        refFaceType
    end

    properties ( SetAccess = protected )
        % number of interpolation nodes associated with the nodal basis function 
        % :math:`\varphi_k(r)`
        nNodal
        % coordinates value of interpolation points
        refNR
        % coordinates value of interpolation points
        refNS
        % coordinates value of interpolation points
        refNT
        % node index of facial points
        FToN
        % number of facial interpolation points
        nFaceNodal
        % total number of face points
        nFaceNodalTotal
    end

    properties ( SetAccess = protected)
        % Vandermonde matrix
        % project from model basis value to nodel basis value
        V
        % % project matrx from interpolation points to quadrature points
        % Vq
        % mass matrix with
        M
        % inverse of mass matrix
        invM
        % derivative matrix with 
        % :math:`[\mathrm{Dr}]_{ij} = 
        % \left.\frac{\partial l_j}{\partial r}\right|_{r_i}`.
        Dr
        Ds
        Dt
        %> lift matrix, \f$ LIFT = M^{-1} \cdot M_e \f$
        LIFT
    end

    methods(Abstract, Access=protected)
        %> get the total number and coordinate of interpolation points
        [ Np,r,s,t ] = getRefNodal(obj, N)
        %> get the derivative of orthogonal function at each interpolation points
        [ dr, ds, dt ] = derivative_orthogonal_func(obj, N, ind, r, s, t);
    end

    methods
        function obj = std_element(N)
            %STD_ELEMENT Construct an instance of this class
            %   Detailed explanation goes here
            obj.num_order = N;
        end

        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.N + inputArg;
        end
    end
end

