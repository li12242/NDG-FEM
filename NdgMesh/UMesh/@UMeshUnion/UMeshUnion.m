%> @brief Super class for unstructed mesh class.
%
%> And here we can put some more detailed informations about the class.
% ======================================================================
%> This class is part of the NDG-FEM software.
%> @author li12242, Tianjin University, li12242@tju.edu.cn
% ======================================================================
classdef UMeshUnion < handle
    
    properties
        %> pointer to cell object of StdCell type
        cell
        %> unstructured edge object
        uedge
        %> number of elements
        K
        %> number of vertices
        Nv
        %> vertices index in each element
        EToV
        %> region index of each element
        EToR @int8
        %> the vertices list of each boundary edge
        BCToV
        %> x component of the vertices coordinates
        vx
        %> y component of the vertices coordinates
        vy
        %> z component of the vertices coordinates
        vz
        %> determination of Jacobian matrix on each IPPS
        J
    end
    
    properties(SetAccess=protected)
%         %> global index of faces adjacent to each element
%         EToFG
%         %> index of element adjacent to each element
%         EToE
%         %> local index of faces adjacent to each element
%         EToF
        %> x coordinates of each IPPS in element
        x
        %> y coordinates of each IPPS in element
        y
        %> z coordinates of each IPPS in element
        z
        %> transformation coefficient of \f$ \partial r/\partial x \f$
        rx
        %> transformation coefficient of \f$ \partial r/\partial y \f$
        ry
        %> transformation coefficient of \f$ \partial r/\partial z \f$
        rz
        %> transformation coefficient of \f$ \partial s/\partial x \f$
        sx
        %> transformation coefficient of \f$ \partial s/\partial y \f$
        sy
        %> transformation coefficient of \f$ \partial s/\partial z \f$
        sz
        %> transformation coefficient of \f$ \partial t/\partial x \f$
        tx
        %> transformation coefficient of \f$ \partial t/\partial y \f$
        ty
        %> transformation coefficient of \f$ \partial t/\partial z \f$
        tz
    end
    
%     properties(SetAccess=protected)
%         %> index of IPPS on edges of each element
%         eidM
%         %> index of adjacent IPPS of each element
%         eidP
%         %> edge type of each IPPS on edges
%         eidtype @int8
%     end
    
%     properties(SetAccess=protected)
%         %> x component of the outward normal vector
%         nx
%         %> y component of the outward normal vector
%         ny
%         %> z component of the outward normal vector
%         nz
%         %> determination of Jacobian matrix on each edge IPPS
%         Js
%     end
    
    properties(Hidden=true)
        %> figure handles
        draw_h
    end
    
    methods
        function obj = UMeshUnion(cell, Nv, vx, vy, vz, K, EToV, EToR, BCToV)
            % set the properties
            obj.cell = cell;
            obj.Nv = Nv;
            obj.vx = vx;
            obj.vy = vy;
            obj.vz = vz;
            obj.K  = K;
            obj.EToV = EToV;
            obj.EToR = int8(EToR);
            obj.BCToV = BCToV;
            
%             [ obj.EToFG ] = getElementalAdjacentFaceGlobalIndex(obj);
%             [ obj.EToE, obj.EToF ] = assembleEleConnect( obj, obj.EToFG );
            [ obj.x, obj.y, obj.z ] = assembleNodeCoordinate(obj, vx, vy, vz);
            
            [obj.rx, obj.ry, obj.rz, obj.sx, obj.sy, obj.sz, ...
                obj.tx, obj.ty, obj.tz, obj.J] ...
                = getElementalNodeInfo(obj);
            
%             [obj.nx, obj.ny, obj.nz, obj.Js] ...
%                 = getElementalSurfaceInfo(obj, obj.vx, obj.vy, obj.vz, obj.EToV);
            
%             [ obj.eidM, obj.eidP, obj.eidtype ] ...
%                 = assembleNodeConnection(obj);
            
            %[ obj.vol, obj.elen ] = getElementalScale(obj);
        end% func
        
        function c_m = getCellIntegralAverage(obj, f_Q)
            % calculate the cell averaged values
            [ c_m ] = cell_mean( f_Q, obj.cell.w, obj.J );
        end
        
        function f_m = getCellFacialIntegralAverage(obj, f_Q)
            % 计算单元内每个面均值
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
        end% func
    end
    
    methods(Abstract, Hidden, Access = protected)
        %> get volume infomation of all the IPPS in each elements
        [ rx, ry, rz, sx, sy, sz, tx, ty, tz, J ] = getElementalNodeInfo(obj)
        %> get out normal vector and surface integral Jacobian coefficient
        %[nx, ny, nz, Js] = getElementalSurfaceInfo(obj)
        %> get the EToFG of each element.
        %[ EToFG ] = getElementalAdjacentFaceGlobalIndex( obj )
    end
    
    methods(Abstract)
        %> refine each element and return a finer mesh object
        refine(obj, order);
        %> draw the mesh discretization
        draw(obj, varargin);
    end
    
    methods(Hidden, Access=protected)
        
%         %> find the adjacent element index and its face index
%         function [EToE, EToF] = assembleEleConnect(obj, EToFG)
%             EToE = ones(obj.cell.Nface, 1)*(1:obj.K);
%             EToF = (1:obj.cell.Nface)'*ones(1,obj.K);
%             
%             for n = 1:(obj.cell.Nface*obj.K)
%                 m = find( abs(EToFG - EToFG(n))<1e-10 );
%                 t = m( m ~= n );
%                 if( ~isempty(t) )
%                     EToE(n) = fix( (t-1)./obj.cell.Nface )+1;
%                     EToF(n) = rem(t-1, obj.cell.Nface)+1;
%                 end
%             end
%         end
        
%         %> calculate the elemental length and the edge length
%         function [ len, eglen ] = getElementalScale(obj)
%             len = sum( bsxfun(@times, obj.cell.w, obj.J) );
%             par = bsxfun(@times, obj.cell.ws, obj.Js);
%             eglen = zeros(obj.cell.Nface, obj.K);
%             if (obj.cell.Nfp == 1)
%                 eglen = par;
%             else
%                 st = 1;
%                 for f = 1:obj.cell.Nface
%                     sk = st + obj.cell.Nfp(f);
%                     row = st:(sk-1);
%                     eglen(f, :) = sum(par(row, :));
%                     st = sk;
%                 end
%             end
%         end% func
        
        %> get the IPPS coordinates in each element
        function [x, y, z] = assembleNodeCoordinate(obj, vx, vy, vz)
            x = obj.proj_vert2node(vx);
            y = obj.proj_vert2node(vy);
            z = obj.proj_vert2node(vz);
        end% func
        
        %> project scalars from verts to nodes in each element
        function nodeQ = proj_vert2node(obj, vertQ)
            ele_vQ = vertQ(obj.EToV);
            nodeQ = obj.cell(1).project_vert2node(ele_vQ);
        end
        
%         function [ idm, idp, type ] = assembleNodeConnection(obj)
%             
%             totalNfp = obj.cell.TNfp;
%             idm = zeros(totalNfp, obj.K);
%             idp = zeros(totalNfp, obj.K);
%             type = int8(zeros(totalNfp, obj.K));
%             
%             Np = obj.cell.Np;
%             for k1 = 1:obj.K
%                 idm(:, k1) = (k1-1)*Np + obj.cell.Fmask(:);
%                 f1 = 1; nfp1 = 1; % 单元内面节点所在面编号与局部节点编号
%                 for n = 1:totalNfp
%                     if nfp1 > obj.cell.Nfp(f1) % 若面节点编号超过f1面节点总数
%                         nfp1 = 1;
%                         f1 = f1 + 1; % 面编号+1
%                     end
%                     type(n, k1) = int8(obj.EToBS(f1, k1));
%                     
%                     k2 = obj.EToE(f1, k1); % 相邻单元
%                     f2 = obj.EToF(f1, k1); % 相邻面
%                     
%                     xpM = obj.x(idm(n, k1)); % 本单元边界节点坐标
%                     ypM = obj.y(idm(n, k1));
%                     zpM = obj.z(idm(n, k1));
%                     % 相邻面上所有节点坐标
%                     ind2 = (k2-1)*Np + obj.cell.Fmask(:, f2);
%                     xP = obj.x( ind2 );
%                     yP = obj.y( ind2 );
%                     zP = obj.z( ind2 );
%                     d12 = (xpM - xP).^2 + (ypM - yP).^2 + (zpM - zP).^2;
%                     m = (d12 < 3e-16);
%                     idp(n, k1) = (k2-1)*Np + obj.cell.Fmask(m, f2);
%                     nfp1 = nfp1 + 1;
%                 end
%             end
%         end
    end
        
    
    
end% classdef

