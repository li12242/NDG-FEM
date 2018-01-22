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
        function obj = StdSubTri(N)
            obj = obj@StdTri(N);
            [ obj.NFV, obj.EToV ] = assembleLocalFiniteVolume( obj );
            [ obj.NCV, obj.LAVToCV ] = assembleLocalControlVolume( obj );
            [ obj.Nedge, obj.n1, obj.n2 ] = assembleLocalEdgeInfo( obj );
            
            [ obj.P, obj.R ] = assembleProjectReconstructMatrix( obj );
            [ obj.VLIFT ] = assembleVLIFT( obj );
            [ obj.SLIFT ] = assembleSLIFT( obj );
        end
        
        function [ nx, ny, nz, Js ] = assembleNormalVector( obj, x, y, z )
            [ nx, ny, nz, Js ] = assembleNormalVector@StdTri( obj, x, y, z );
            Js = zeros( size(Js) );
            
            for i1 = 1:(obj.TNfp-1)
                i2 = i1+1;
                v1 = obj.Fmask(i1);
                v2 = obj.Fmask(i2);
                d2 = sqrt( (x(v1, :) - x(v2, :)).^2 ...
                    + (y(v1, :) - y(v2, :)).^2 )/2;
                Js(i1, :) = Js(i1, :) + d2;
                Js(i2, :) = Js(i2, :) + d2;
            end
        end
        
        function [nxf, nyf, nzf, flen] = ...
                matEvaluateInnerControlEdgeInfo( obj, x, y, z )
            K = size(x, 2);
            nxf = zeros( obj.Nedge, K );
            nyf = zeros( obj.Nedge, K );
            nzf = zeros( obj.Nedge, K );
            flen = zeros( obj.Nedge, K );
            
            xc = zeros( obj.Nedge, K );
            yc = zeros( obj.Nedge, K );
            for k = 1:K
                temp = x(:, k);
                temp1 = repmat( mean( temp(obj.EToV) ), obj.Nface, 1 );
                xc(:,k) = temp1(:);
                temp = y(:, k);
                temp1 = repmat( mean( temp(obj.EToV) ), obj.Nface, 1 );
                yc(:,k) = temp1(:);
            end
            for k = 1:K
                x1 = x( obj.n1, k ); 
                y1 = y( obj.n1, k );
                x2 = x( obj.n2, k ); 
                y2 = y( obj.n2, k );
                xf = 0.5*( x1 + x2 ); 
                yf = 0.5*( y1 + y2 );
                nxt = yc - yf; 
                nyt = xf - xc;
                ds = sqrt( nxt.^2 + nyt.^2 );
                nxf( :, k ) = nxt./ds;
                nyf( :, k ) = nyt./ds;
                flen( :, k ) = ds;
            end
        end% func
        
        function draw(obj, n)
            
        end
    end
    
    methods( Hidden, Access = protected )
        %> Get the total number of FV and the node index in each FV
        function [ NFV, EToV ] = assembleLocalFiniteVolume( obj )
            NFV = obj.N^2;
            EToV = zeros(3, NFV);

            sk = 1;
            for row = 1:obj.N
                v1 = obj.Fmask(row, 3);
                v2 = obj.Fmask(row+1, 3);
                for kb = 1:row
                    EToV(:, sk) = [v1, v2, v2+1]';
                    v1 = v1+1;
                    v2 = v2+1;
                    sk = sk+1;
                end
                
                v1 = obj.Fmask(row, 3);
                v2 = obj.Fmask(row+1, 3);
                for kt = 1:row-1
                    EToV(:, sk) = [v1, v2+1, v1+1]';
                    v1 = v1+1;
                    v2 = v2+1;
                    sk = sk+1;
                end
            end
        end
        
        %> Get the total number of CV and the length/area/volume of each CV
        function [ NCV, LAVToCV ] = assembleLocalControlVolume( obj )
            NCV = obj.Np;
            LAVToCV = zeros(NCV, 1);
            for k = 1:obj.NFV
                vert = obj.EToV(:, k); % get the node index in each FV
                area = TriArea( obj.r(vert), obj.s(vert) );
                LAVToCV(vert) = LAVToCV(vert) + area/3;
            end
        end% func
        
        %> Get the total number of edge, the outward normal vector and the
        %> adjacent nodes index for each edge.
        function [ Nedge, v1, v2 ] = assembleLocalEdgeInfo( obj )
            Nface = obj.Nface;
            Nedge = obj.NFV*Nface;
            v1 = zeros(Nface, obj.NFV);
            v2 = zeros(Nface, obj.NFV);
            for k = 1:obj.NFV % loop over FVs
                vertId = obj.EToV(:, k);
                for f1 = 1:Nface
                    f2 = mod(f1,Nface)+1;
                    fv1 = vertId( f1 );
                    fv2 = vertId( f2 );
                    
                    v1(f1, k) = fv1;
                    v2(f1, k) = fv2;
                end
            end
            v1 = v1(:);
            v2 = v2(:);
        end
        
        %> Get the project and reconstruction matrix
        [ P, R ] = assembleProjectReconstructMatrix( obj )
        
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
        
    end% methods
    
end
