classdef StdSubTri < StdTri
    
    properties
        %> No. of finite volume
        NFV
        %> vertex index in each finite volume
        EToV
        %> No. of control volume, equal to Np
        NCV
        %> length/area/volume of each control volume
        LAVToCV
        %> No. of inside edges, equal to NFV*3
        Nedge
        %> outward normal vector of each edge, from v1 to v2
        nr, ns
        %> Adjacent nodes index of each edge
        n1, n2
        VLIFT
        %> 
        SLIFT
        %> project matrix to project basis parameters to coutrol volume values
        P
        %> reconstruct matrix, satisfying P*R = I
        R
    end
    
    methods
        function obj = StdSubTri(N)
            obj = obj@StdTri(N);
            [ obj.NFV, obj.EToV ] = assembleLocalFiniteVolume( obj );
            [ obj.NCV, obj.LAVToCV ] = assembleLocalControlVolume( obj );
            [ obj.Nedge, obj.nr, obj.ns, obj.n1, obj.n2 ] = assembleLocalEdge( obj );
            obj.VLIFT = zeros( obj.Np, obj.Nedge );
            for e = 1:obj.Nedge
                n1 = obj.n1(e);
                n2 = obj.n2(e);
                obj.VLIFT(n1, e) = -1;
                obj.VLIFT(n2, e) =  1;
            end
            [ obj.P, obj.R ] = assembleProjectReconstructMatrix( obj );
            [ obj.SLIFT ] = assembleLIFT( obj );
        end
        
        
    end
    
    methods( Hidden, Access = protected )
        %> Get the total number of FV and the node index in each FV
        [ NFV, EToV ] = assembleLocalFiniteVolume( obj )
        %> Get the total number of CV and the length/area/volume of each CV
        [ NCV, LAV ] = assembleLocalControlVolume( obj )
        %> Get the total number of edge, the outward normal vector and the
        %> adjacent nodes index for each edge.
        [ Nedge, nx, ny, v1, v2 ] = assembleLocalEdge( obj )
        %> Get the project and reconstruction matrix
        [ P, R ] = assembleProjectReconstructMatrix( obj )
        
        function [ LIFT ] = assembleLIFT(obj)
            LIFT = zeros(obj.Np, obj.TNfp);
            for n = 1:obj.TNfp
                row = obj.Fmask(n);
                LIFT(row, n) = 1.0;
            end
        end
    end
    
end

