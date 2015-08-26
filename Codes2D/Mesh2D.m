classdef Mesh2D < handle
    properties
        Element     % element object
        
        K           % element number
        VX          % vertices x coordinate
        VY          % vertices y coordinate
        Fx          % node point coordinate
        Fy          % face node coordinate: x
        x           % face node coordinate: y
        y
        rx
        ry
        sx
        sy
        J
        nx
        ny
        sJ
        EToE
        EToF
%         mapM 
        vmapM       % face node to Inner node
        vmapP       % face node to Outer node
        mapP        % inverse of vmapP
        vmapB       % face node to Boundary node
        mapB        % boundary node to face node
        BCType
        mapI,   vmapI
        mapO,   vmapO
        mapW,   vmapW
        mapF,   vmapF
        mapC,   vmapC
        mapN,   vmapN
        mapD,   vmapD
        mapS,   vmapS
        
        % Low storage Runge-Kutta coefficients
        rk4a = [            0.0 ...
                -567301805773.0/1357537059087.0 ...
                -2404267990393.0/2016746695238.0 ...
                -3550918686646.0/2091501179385.0  ...
                -1275806237668.0/842570457699.0];
        rk4b = [ 1432997174477.0/9575080441755.0 ...
                 5161836677717.0/13612068292357.0 ...
                 1720146321549.0/2090206949498.0  ...
                 3134564353537.0/4481467310338.0  ...
                 2277821191437.0/14882151754819.0];
        rk4c = [             0.0  ...
                 1432997174477.0/9575080441755.0 ...
                 2526269341429.0/6820363962896.0 ...
                 2006345519317.0/3224310063776.0 ...
                 2802321613138.0/2924317926251.0]; 
       
    end %properties
    
    methods
        function obj = Mesh2D(Element, EToV, VX, VY, BCType)
            obj.K = size(EToV, 1);
            obj.VX = VX; obj.VY = VY;
            obj.Element = Element;
            va = EToV(:,1)'; vb = EToV(:,2)'; vc = EToV(:,3)';
            obj.x = 0.5*(-(Element.r+Element.s)*VX(va)+(1+Element.r)*VX(vb)+(1+Element.s)*VX(vc));
            obj.y = 0.5*(-(Element.r+Element.s)*VY(va)+(1+Element.r)*VY(vb)+(1+Element.s)*VY(vc));
            obj.Fx = obj.x(Element.Fmask(:), :); obj.Fy = obj.y(Element.Fmask(:), :);
            [obj.rx,obj.sx,obj.ry,obj.sy,obj.J] = GeometricFactors2D(obj.x,obj.y,Element.Dr,Element.Ds);
            [obj.nx, obj.ny, obj.sJ] = Normals2D(Element, obj);
            [obj.EToE, obj.EToF] = tiConnect2D(EToV);
            BuildMaps2D(Element,obj,EToV);
            obj.BCType = BCType;
        end %function
    end %methods
    
end %class