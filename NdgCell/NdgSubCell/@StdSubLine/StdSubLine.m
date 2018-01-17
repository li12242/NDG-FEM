classdef StdSubLine < StdLine
    %STDSUBLINE Summary of this class goes here
    %   Detailed explanation goes here
    
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
        %> local nodes index of each edge
        n1
        %> adjacent nodes index of each edge
        n2
        %> volume lift matrix to add inner edge fluxes to RHS terms
        VLIFT
        %> surface lift matrix to add surface edge fluxes to RHS terms
        SLIFT
        %> project matrix to project basis parameters to coutrol volume values
        P
        %> reconstruct matrix, satisfying P*R = I
        R
    end
    
    methods
        
        function obj = StdSubLine(N)
            obj = obj@StdLine(N);
            [ obj.NFV, obj.EToV ] = assembleLocalFiniteVolume( obj );
            [ obj.NCV, obj.LAVToCV ] = assembleLocalControlVolume( obj );
            [ obj.Nedge, obj.n1, obj.n2 ] = assembleLocalEdgeInfo( obj );
            
            [ obj.P, obj.R ] = assembleProjectReconstructMatrix( obj );
            [ obj.VLIFT ] = assembleVLIFT( obj );
            [ obj.SLIFT ] = assembleSLIFT( obj );
        end
        
        function [ nx, ny, nz, Js ] = assembleNormalVector( obj, x, y, z )
            [ nx, ny, nz, Js ] = assembleNormalVector@StdLine( obj, x, y, z );
            Js = ones( size(Js) );            
        end
        
        function [nxf, nyf, nzf, flen] = ...
                matEvaluateInnerControlEdgeInfo( obj, x, y, z )
            K = size(x, 2);
            
            nxf = zeros( obj.Nedge, K );
            nyf = zeros( obj.Nedge, K );
            nzf = zeros( obj.Nedge, K );
            flen = ones( obj.Nedge, K );
            
            for k = 1:K
                x1 = x( obj.n1, k ); 
                x2 = x( obj.n2, k ); 
                nxf( :, k ) = sign( x2 - x1 );
            end
        end
        
        function draw(obj, n)
            
        end
    end
    
    methods( Hidden, Access = protected )
        %> Get the total number of FV and the node index in each FV
        function [ NFV, EToV ] = assembleLocalFiniteVolume( obj )
            NFV = obj.N;
            EToV = zeros(2, NFV);
            EToV(1, :) = (1:obj.N);
            EToV(2, :) = EToV(1, :) + 1;
        end
        
        %> Get the total number of CV and the length/area/volume of each CV
        function [ NCV, LAVToCV ] = assembleLocalControlVolume( obj )
            NCV = obj.Np;
            LAVToCV = zeros(NCV, 1);
            for k = 1:obj.NFV
                vert = obj.EToV(:, k);
                len = abs( diff( obj.r( vert ) ) );
                LAVToCV(vert) = LAVToCV(vert) + len/2;
            end
            %LAVToCV = LAVToCV./obj.LAV;
        end
        
        %> Get the total number of edge, the outward normal vector and the
        %> adjacent nodes index for each edge.
        function [ Nedge, v1, v2 ] = assembleLocalEdgeInfo( obj )
            Nedge = obj.N;
            v1 = obj.EToV(1, :);
            v2 = obj.EToV(2, :);
            v1 = v1(:); 
            v2 = v2(:);
        end
        
        %> Get the project and reconstruction matrix
        function [ P, R ] = assembleProjectReconstructMatrix( obj )
            Nface = obj.Nface;
            r = obj.r;
            P = zeros( obj.Np, obj.Np );
            Pt = zeros( obj.Nq, obj.Np);
            for k = 1:obj.NFV
                vertId = obj.EToV(:, k);
                rc = mean( r(vertId) );
                
                for f = 1:Nface
                    v1 = vertId(f);
                    r1 = r(v1);
                    rq1 = obj.project_vert2quad([r1, rc]');
                    fvArea = abs( diff( [r1, rc] ) );
                    for i = 1:obj.Np
                        Pt(:,i) = obj.orthogonal_func(obj.N, i, rq1);
                    end
                    P(v1, :) = P(v1, :) + (fvArea .* obj.wq' ./ obj.LAV)*Pt;
                end
            end
            P = P/obj.V;
            P = bsxfun(@times, P, 1./obj.LAVToCV);
            R = inv(P);
        end
        
        function [ VLIFT ] = assembleVLIFT( obj )
            VLIFT = zeros( obj.Np, obj.Nedge );
            ind = sub2ind( [obj.Np, obj.Nedge], obj.n1', 1:obj.Nedge );
            VLIFT( ind ) = -1;
            ind = sub2ind( [obj.Np, obj.Nedge], obj.n2', 1:obj.Nedge );
            VLIFT( ind ) = 1;
        end
        
        function [ SLIFT ] = assembleSLIFT(obj)
            SLIFT = zeros(obj.Np, obj.TNfp);
            for n = 1:obj.TNfp
                row = obj.Fmask(n);
                SLIFT(row, n) = 1.0;
            end
        end% func
    end
    
end

