classdef Element1D < handle
    properties
        N       % degree
        Np      % node point
        Nfp = 1     % number of point at face
        Nfaces = 2  % face number
        r       % local coordinate
        V       %
        invV
        Dr
        Fmask   % edge nodes to element nodes
        LIFT    % M^{-1}*M_{edge}
        
    end %properties
    
    methods
        function obj = Element1D(N)
            NODETOL = 10e-5;
            
            obj.N = N; obj.Np = N+1; obj.Nfp = 1; obj.Nfaces=2;
            obj.r = JacobiGL(0,0,N);
            obj.V  = Vandermonde1D(N, obj.r); obj.invV = inv(obj.V);
            obj.Dr = Dmatrix1D(N, obj.r, obj.V);
            obj.LIFT = Lift1D(obj);
            
            fmask1 = find( abs(obj.r+1) < NODETOL)';
            fmask2 = find( abs(obj.r-1) < NODETOL)';
            obj.Fmask  = [fmask1;fmask2]';
            
        end
    end
end %classdef
    