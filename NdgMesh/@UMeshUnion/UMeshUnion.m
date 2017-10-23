%> @brief Here we have a brief description of the class.
%
%> And here we can put some more detailed informations about the class.
% ======================================================================
%> This class is part of the NDGOM software.
%> @author li12242, Tianjin University, li12242@tju.edu.cn
% ======================================================================
classdef UMeshUnion < handle
    properties
        %> standard cell class
        cell
        %> number of cell
        K
        %> number of vertices
        Nv
        %> vertex index in each cell (column)
        EToV
        %> region types for each cell
        EToR
        %> boundary types for each cell (column)
        EToB
        %> coordinate of vertex
        vx
        %> coordinate of vertex
        vy
        %> coordinate of vertex
        vz
    end
    
    % elemental volume infomation
    properties( Hidden, SetAccess=protected )
        %> mesh id of adjacent cell
        EToM
        %> adjacent cell index for each cell
        EToE
        %> adjacent face index for each cell
        EToF
        %> coordinate of interpolation points
        x
        %> coordinate of interpolation points
        y
        %> coordinate of interpolation points
        z
        %> determination of Jacobian matrix at each interpolation points
        J
        rx, ry, rz
        sx, sy, sz
        tx, ty, tz
        %> length/area/volume of each cell
        vol
        %> character length of each cell
        elen
        %spg_delta
    end
    
    properties( Hidden, SetAccess = protected )
        %> the local face point index
        eidM
        %> the adjacent point index
        eidP
        %> edge type of each facial points
        eidtype
    end
    
    % elemental edge information
    properties( Hidden, SetAccess = protected )
        %> normal vector of each facial point
        nx
        %> normal vector of each facial point
        ny
        %> normal vector of each facial point
        nz
        %> determination of facial integral at each face points
        Js
    end
    
    properties( Hidden = true )
        figureHandle  % figure handle
    end
    
    methods( Abstract, Hidden, Access = protected )
        assembleJacobiFactor(obj) % get volume infomation () of each element
        assembleFacialJaobiFactor(obj) % get out normal vector of each elemental edges
        faceId = assembleGlobalFaceIndex(obj)
    end
    
    methods(Hidden, Access=protected)
        [ EToE, EToF, EToM ] = assembleCellConnect(obj)
        [ eidM, eidP, eidtype ] = assembleEdgeNode(obj)
        [ EToB ] = assembleCellBoundary( obj, BCToV )
        
        function obj = ele_scale(obj)
            len = sum( bsxfun(@times, obj.cell.w, obj.J) );
            par = bsxfun(@times, obj.cell.ws, obj.Js);
            eglen = zeros(obj.cell.Nface, obj.K);
            if (obj.cell.Nfp == 1)
                eglen = par;
            else
                st = 1;
                for f = 1:obj.cell.Nface
                    sk = st + obj.cell.Nfp(f);
                    row = st:(sk-1);
                    eglen(f, :) = sum(par(row, :));
                    st = sk;
                end
            end
            obj.vol = len;
            obj.elen = eglen;
        end% func
        
        function [x, y, z] = assembleNodeCoor(obj, vx, vy, vz)
            x = obj.proj_vert2node(vx);
            y = obj.proj_vert2node(vy);
            z = obj.proj_vert2node(vz);
        end% func
    end
    
    methods
        
        function obj = UMeshUnion(cell, Nv, vx, vy, vz, K, EToV, EToR, BCToV)
            [obj.cell, obj.Nv, obj.vx, obj.vy, obj.vz, obj.K, obj.EToV, obj.EToR] ...
                = checkInput(cell, Nv, vx, vy, vz, K, EToV, EToR);
                        
            [ obj.EToE, obj.EToF, obj.EToM]  = assembleCellConnect( obj );
            [ obj.x, obj.y, obj.z ] = assembleNodeCoor( obj, vx, vy, vz );
            [ obj.EToB ] = assembleCellBoundary(obj, BCToV);
            
            [obj.rx, obj.ry, obj.rz, obj.sx, obj.sy, obj.sz, ...
                obj.tx, obj.ty, obj.tz, obj.J] = assembleJacobiFactor( obj );
            [obj.nx, obj.ny, obj.nz, obj.Js] ...
                = assembleFacialJaobiFactor( obj );
            
            [obj.eidM, obj.eidP, obj.eidtype] = assembleEdgeNode( obj );
            
            %obj = ele_scale(obj);
        end% func
        
        function nodeQ = proj_vert2node(obj, vertQ)
            % project scalars from mesh verts to nodes
            ele_vQ = vertQ(obj.EToV);
            nodeQ = obj.cell.project_vert2node(ele_vQ);
        end
        
        function c_m = cell_mean(obj, f_Q)
            % calculate the cell averaged values
            [ c_m ] = cell_mean( f_Q, obj.cell.w, obj.J );
        end
        
        function f_m = face_mean(obj, f_Q)
            par = bsxfun(@times, obj.cell.ws, obj.Js);
            f_m = zeros(obj.cell.Nface, obj.K);
            f_M = f_Q(obj.eidM);
            tmp = par.*f_M;
            st = 1;
            if (obj.cell.Nfp == 1)
                f_m = tmp;
            else
                for f = 1:obj.cell.Nface
                    sk = st + obj.cell.Nfp(f);
                    row = st:(sk-1);
                    f_m(f, :) = sum( tmp(row, :) )./sum(par(row, :));
                    st = sk;
                end
            end
        end
    end
    
end

function [cell, Nv, vx, vy, vz, K, EToV, EToR] ...
    = checkInput(cell, Nv, vx, vy, vz, K, EToV, EToR)

if( numel(vx) ~= Nv ) || (numel(vy) ~= Nv ) || (numel(vz) ~= Nv )
    msgID = 'UMeshUnion:InputVertexError';
    msgtext = 'The number of input vertex is not equal to Nv.';
    ME = MException(msgID, msgtext);
    throw(ME);
end

if isrow(vx) 
    vx = vx'; 
end

if isrow(vy) 
    vy = vy'; 
end

if isrow(vz) 
    vz = vz'; 
end

if( size(EToV, 1) ~= cell.Nv )
    msgID = 'UMeshUnion:InputCellError';
    msgtext = 'The rows of EToV is not equal to Nv (cell).';
    ME = MException(msgID, msgtext);
    throw(ME);
end

if( size(EToV, 2) ~= K ) || (numel(EToR) ~= K)
    msgID = 'UMeshUnion:InputCellError';
    msgtext = 'The number of cells in EToV or EToR is not equal to K.';
    ME = MException(msgID, msgtext);
    throw(ME);
end

EToR = NdgRegionType( EToR );
end

