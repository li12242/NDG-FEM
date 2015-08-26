classdef Mesh1D < handle
    properties
        Element   % element object
        
        K           % element number
        VX          % vertices coordinate
        x           % node coordinate
        Fx          % vertices coordinate (What is the different with VX?)
        nx          % outward pointing normals
        
        Fscale      % jacobi_{edge}/jacobi_{volume}
        J           % jacobi_{volume}
        rx          % $\partial r/ \partial x$
        vmapM       % face node to Inner node
        vmapP       % face node to Outer node 
        vmapB       % face node to Boundary node
        mapB        % boundary node to face node
        vmapI       % inner face node to node
        vmapO       % outer face node to node
        mapI        % inner boundary node to face node
        mapO        % outer boundary node to face node
        
        
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
    end %properties public
    
    properties(SetAccess=private, Hidden, GetAccess=private)
        EToE        % Element to Element
        EToF        % Element to Face
    end %properties provate
    
    methods
        function obj = Mesh1D(Element, EToV, VX)
            obj.Element = Element;
            obj.VX = VX;
            
            va = EToV(:,1)'; vb = EToV(:,2)';
            obj.K = size(EToV, 1);
            obj.x = ones(Element.N+1,1)*VX(va) + 0.5*(Element.r+1)*(VX(vb)-VX(va));
            
            [obj.rx,obj.J] = GeometricFactors1D(obj.x, Element.Dr);
            obj.Fx = obj.x(Element.Fmask(:), :);
            
            % Build surface normals and inverse metric at surface
            [obj.nx] = Normals1D(obj);
            obj.Fscale = 1./(obj.J(Element.Fmask,:));
            [obj.EToE, obj.EToF] = Connect1D(EToV);
            
            obj = BuildMaps1D(obj);
        end
        
        function Mesh1D = BuildMaps1D(Mesh1D)

            % function [vmapM, vmapP, vmapB, mapB] = BuildMaps1D
            % Purpose: Connectivity and boundary tables for nodes given in the K # of elements,
            % 	       each with N+1 degrees of freedom.

            % Globals1D;
            NODETOL = 10e-5;

            K = Mesh1D.K;  Nfp = Mesh1D.Element.Nfp;
            Np = Mesh1D.Element.Np; Nfaces = Mesh1D.Element.Nfaces;

            % number volume nodes consecutively
            nodeids = reshape(1:K*Np, Np, K);
            vmapM   = zeros(Nfp, Nfaces, K);
            vmapP   = zeros(Nfp, Nfaces, K);

            for k1=1:K
              for f1=1:Nfaces
                % find index of face nodes with respect to volume node ordering
                vmapM(:,f1,k1) = nodeids(Mesh1D.Element.Fmask(:,f1), k1);
              end
            end

            for k1=1:K
              for f1=1:Nfaces
                % find neighbor
                k2 = Mesh1D.EToE(k1,f1); f2 = Mesh1D.EToF(k1,f1);

                % find volume node numbers of left and right nodes 
                vidM = vmapM(:,f1,k1); vidP = vmapM(:,f2,k2);

                x1  = Mesh1D.x(vidM); x2  = Mesh1D.x(vidP);

                % Compute distance matrix
                D = (x1 -x2 ).^2;
                if (D<NODETOL) vmapP(:,f1,k1) = vidP; end;
              end
            end

            Mesh1D.vmapP = vmapP(:); Mesh1D.vmapM = vmapM(:);

            % Create list of boundary nodes
            Mesh1D.mapB = find(vmapP==vmapM); Mesh1D.vmapB = vmapM(Mesh1D.mapB);

            % Create specific left (inflow) and right (outflow) maps
            Mesh1D.mapI = 1; Mesh1D.mapO = K*Nfaces; Mesh1D.vmapI = 1; Mesh1D.vmapO = K*Np;
            return
        end %function BuildMaps1D
    end
    
end

