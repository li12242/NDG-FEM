%> Create the 2d mesh object. The object is inherited from the mesh
%> object and has special treatment for the 2d elements (triangle and
%> quadrilateral). The mesh object contains only one type of elements,
%> and allocate multiple choices to create the object. The input
%> methods include
classdef NdgMesh2d < NdgMesh
    
    properties( Hidden=true )
        point_h
    end
    
    methods(Hidden, Access=protected)
        function [ edge ] = makeConnectNdgEdge( obj, mesh1, mid0, mid1 )
            edge = NdgEdge2d( obj, mesh1, mid0, mid1 );
        end
        
        function [rx, ry, rz, sx, sy, sz, tx, ty, tz, J] ...
                = assembleJacobiFactor(obj)
            xr = obj.cell.Dr*obj.x; xs = obj.cell.Ds*obj.x;
            yr = obj.cell.Dr*obj.y; ys = obj.cell.Ds*obj.y;
            J = -xs.*yr + xr.*ys;
            
            rx = ys./J; sx =-yr./J;
            ry =-xs./J; sy = xr./J;
            
            tx = zeros(size(rx));
            ty = zeros(size(ry));
            rz = zeros(size(rx));
            sz = zeros(size(rx));
            tz = zeros(size(rx));
        end
        
        function [nx, ny, nz, Js] ...
                = assembleFacialJaobiFactor( obj )
            Nface = obj.cell.Nface;
            TNfp = obj.cell.TNfp;
            nx = zeros(TNfp, obj.K);
            ny = zeros(TNfp, obj.K);
            nz = zeros(TNfp, obj.K);
            
            faceIndexStart = ones(obj.cell.Nface, 1); % start index of each face node
            for f = 2:obj.cell.Nface
                faceIndexStart(f) = faceIndexStart(f-1) + obj.cell.Nfp(f-1);
            end
            
            for f = 1:Nface
                Nfp = obj.cell.Nfp(f);
                face_x1 = obj.vx( obj.EToV( obj.cell.FToV(1,f), :))';
                face_x2 = obj.vx( obj.EToV( obj.cell.FToV(2,f), :))';
                face_y1 = obj.vy( obj.EToV( obj.cell.FToV(1,f), :))';
                face_y2 = obj.vy( obj.EToV( obj.cell.FToV(2,f), :))';
                
                ind = faceIndexStart(f):(faceIndexStart(f)+Nfp-1);
                nx(ind, :) = repmat( (face_y2 - face_y1), Nfp, 1 );
                ny(ind, :) = repmat(-(face_x2 - face_x1), Nfp, 1 );
            end
            Js = sqrt(nx.*nx+ny.*ny);
            % normalise
            nx = nx./Js;
            ny = ny./Js;
            Js = Js.*0.5;
        end
        
        function faceId = assembleGlobalFaceIndex(obj)
            faceId = zeros(obj.cell.Nface, obj.K);
            for f = 1:obj.cell.Nface
                v1 = obj.EToV(obj.cell.FToV(1,f), :);
                v2 = obj.EToV(obj.cell.FToV(2,f), :);
                % calculate the indicator for each edge
                faceId(f, :) = min(v1, v2)*obj.Nv + max(v1, v2);
            end
        end
    end% methods
    
    % public methods
    methods(Access=public)
        function obj = NdgMesh2d(cell, Nv, vx, vy, K, EToV, EToR, BCToV)
            if (nargin ~= 8)
                msgID = [mfilename, ':InputError'];
                msgtext = 'The number of inputs should be 8.';
                ME = MException(msgID, msgtext);
                throw(ME);
            end
            [cell, Nv, vx, vy, K, EToV, EToR, BCToV] ...
                = checkInput(cell, Nv, vx, vy, K, EToV, EToR, BCToV);
            [ EToV ] = makeCounterclockwiseVertexOrder( EToV, vx, vy );
            vz = zeros( size(vx) ); % vz is all zeros
            obj = obj@NdgMesh(cell, Nv, vx, vy, vz, K, EToV, EToR, BCToV);
        end% func
        
        %         function obj = add_sponge(obj, vertlist)
        %             % Calculate the distance from the boundary for each nodes
        %             obj. spg_delta = zeros(obj.cell.Np, obj.K);
        %             xb = obj.vx(vertlist);
        %             yb = obj.vy(vertlist);
        %             for k = 1:obj.K
        %                 if obj.EToR(k) ~= ndg_lib.mesh_type.Sponge
        %                     continue;
        %                 end
        %
        %                 for n = 1:obj.cell.Np
        %                     xi = obj.x(n, k);
        %                     yi = obj.y(n, k);
        %                     obj.spg_delta(n, k) = ...
        %                         min( sqrt( (xi - xb).^2 + (yi - yb).^2 ) );
        %                 end
        %             end
        %         end
        
        %         function spg_sigma = cal_sponge_strength(obj, spg_len, max_sigma)
        %             % ���㺣������ɳ�ϵ��?sigma
        %             spg_sigma = zeros(obj.cell.Np, obj.K);
        %             p = 3;
        %             if isempty(obj.spg_delta)
        %                 return;
        %             end
        %             for k = 1:obj.K
        %                 if obj.EToR(k) ~= ndg_lib.mesh_type.Sponge
        %                     continue;
        %                 end
        %
        %                 spg_sigma(:, k) = ...
        %                     max_sigma*(1 - obj.spg_delta(:, k)/spg_len).^p;
        %             end
        %         end
    end% methods
    
end

function [cell, Nv, vx, vy, K, EToV, EToR, BCToV] ...
    = checkInput(cell, Nv, vx, vy, K, EToV, EToR, BCToV)
% check the input variables for initlizing the mesh object.
if( ~isa(cell, 'StdTri') && ~isa(cell, 'StdQuad') )
    msgID = [mfilename, ':InputStdCellError'];
    msgtext = 'The input standard cell should be a StdTri or StdQuad object.';
    ME = MException(msgID, msgtext);
    throw(ME);
end

if( size(EToV, 1) ~= cell.Nv )
    msgID = [mfilename, ':InputCellError'];
    msgtext = 'The rows of EToV is not equal to Nv (cell).';
    ME = MException(msgID, msgtext);
    throw(ME);
end

if( size(BCToV, 1) ~= 3 )
    msgID = [mfilename, ':InputBCToV'];
    msgtext = 'The rows of input BCToV should be 3 ( [v1, v2, bcType] ).';
    ME = MException(msgID, msgtext);
    throw(ME);
end

EToR = NdgRegionType( EToR );
end% func
