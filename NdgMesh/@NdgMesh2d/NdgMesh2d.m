%> Create the 2d mesh object. The object is inherited from the mesh
%> object and has special treatment for the 2d elements (triangle and
%> quadrilateral). The mesh object contains only one type of elements,
%> and allocate multiple choices to create the object. The input
%> methods include
classdef NdgMesh2d < NdgMesh
    properties( Constant )
        type = enumMeshDim.Two
    end
    
    methods( Hidden, Access = protected )
        
        function obj = assembleJacobiFactor(obj)
            Np = obj.cell.Np;
            K = obj.K;
            
            xr = obj.cell.Dr * obj.x;
            xs = obj.cell.Ds * obj.x;
            yr = obj.cell.Dr * obj.y;
            ys = obj.cell.Ds * obj.y;
            J = -xs .* yr + xr .* ys;
            
            obj.rx = ys ./ J; obj.sx = - yr ./ J;
            obj.ry =-xs ./ J; obj.sy =  xr ./ J;
            
            obj.J = J;
            obj.tx = zeros(Np, K);
            obj.ty = zeros(Np, K);
            obj.rz = zeros(Np, K);
            obj.sz = zeros(Np, K);
            obj.tz = zeros(Np, K);
        end
    end% methods
    
    % public methods
    methods(Access=public)
        function obj = NdgMesh2d(cell, Nv, vx, vy, K, EToV, EToR)
            [ EToV ] = makeCounterclockwiseVertexOrder( EToV, vx, vy );
            vz = zeros( Nv, 1 ); % vz is all zeros
            obj = obj@NdgMesh(cell, Nv, vx, vy, vz, K, EToV, EToR);
        end% func
    end% methods
    
end