classdef Triangle < handle
    properties
        N       % degree
        Np      % number of point
        Nfp     % number of point on edge
        Nfaces  % number of faces
        r       % coordinate
        s       % coordinate
        V       % Vandermonde Matrix
        MassMatrix  % Mass Matrix
        Dr
        Ds      % difference Matrix
        Fmask   % face point number to node number
        LIFT    % surface integral
        Drw     % weak difference matrix
        Dsw
    end %properties
    
    methods
        function obj = Triangle(N)
            NODETOL = 1e-12;
            obj.N = N;
            obj.Np = (N+1)*(N+2)/2;
            obj.Nfp = N+1;
            obj.Nfaces = 3;
            [x,y] = Nodes2D(N); [obj.r,obj.s] = xytors(x,y);
            obj.V = Vandermonde2D(N,obj.r,obj.s); invV = inv(obj.V);
            obj.MassMatrix = invV'*invV;
            [obj.Dr,obj.Ds] = Dmatrices2D(N, obj.r, obj.s, obj.V);
            % find all the nodes that lie on each edge
            fmask1   = find( abs(obj.s+1) < NODETOL)'; 
            fmask2   = find( abs(obj.r+obj.s) < NODETOL)';
            fmask3   = find( abs(obj.r+1) < NODETOL)';
            obj.Fmask  = [fmask1;fmask2;fmask3]';
            obj.LIFT = Lift2D(obj);
            [Vr, Vs] = GradVandermonde2D(N, obj.r, obj.s);
            obj.Drw = (obj.V*Vr')/(obj.V*obj.V'); 
            obj.Dsw = (obj.V*Vs')/(obj.V*obj.V');
        end
    end %methods
end %classdef