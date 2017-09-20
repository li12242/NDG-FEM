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
            
            [ obj.x, obj.y, obj.z ] = assembleNodeCoordinate(obj, vx, vy, vz);
            
            [obj.rx, obj.ry, obj.rz,obj.sx, obj.sy, obj.sz, ...
                obj.tx, obj.ty, obj.tz, obj.J] = getElementalNodeInfo(obj);
            
            [ obj.uedge ] = setUEdgeClass(obj, BCToV);
            %[ obj.vol, obj.elen ] = getElementalScale(obj);
        end% func
        
        %> calculate the cell integral averaged values
        function c_m = getCellIntegralAverage(obj, f_Q)
            [ c_m ] = cell_mean( f_Q, obj.cell.w, obj.J );
        end
        
        %> Calculate the facial averages of each element
        function f_m = getCellFacialIntegralAverage(obj, f_Q)
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
        [ uedge ] = setUEdgeClass(obj, BCToV);
    end
    
    methods(Abstract)
        %> refine each element and return a finer mesh object
        refine(obj, order);
        %> draw the mesh discretization
        draw(obj, varargin);
    end
    
    methods(Hidden, Access=protected)
        
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
        
    end
    
    
    
end% classdef

