classdef NdgMesh1d < NdgMesh
    %LINE_MESH Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(Constant)
        type = enumMeshDim.One
    end
    
    methods(Hidden, Access = protected)
        function [ edge ] = makeConnectNdgEdge( obj, mesh1, mid0, mid1 )
            edge = NdgEdge1d( obj, mesh1, mid0, mid1 );
        end
        
        function obj = assembleJacobiFactor( obj )
            xr = obj.cell.Dr * obj.x;
            obj.J = xr; 
            obj.rx = 1./obj.J;
            
%             [r, c] = size(rx);
%             obj.ry = zeros(r, c);
%             obj.rz = zeros(size(rx));
%             obj.sx = zeros(size(rx));
%             obj.sy = zeros(size(rx));
%             obj.sz = zeros(size(rx));
%             obj.tx = zeros(size(rx));
%             obj.ty = zeros(size(rx));
%             obj.tz = zeros(size(rx));
        end
        
%         function [nx, ny, nz, Js] = assembleFacialJaobiFactor( obj )
%             % Define outward normals
%             xb = obj.x(obj.cell.Fmask, :);
%             xc = obj.GetMeshAverageValue( obj.x );
%             nx = sign( xb - xc );
%             
%             Js = ones(size(nx));
%             ny = zeros(obj.cell.Nface, obj.K);
%             nz = zeros(obj.cell.Nface, obj.K);
%         end
        
%         function faceId = assembleGlobalFaceIndex( obj )
%             faceId = zeros(obj.cell.Nface, obj.K);
%             for f = 1:obj.cell.Nface
%                 faceId(f, :) = obj.EToV(obj.cell.FToV(1,f), :);
%             end
%         end% func
    end% methods
    
    methods
        %> \brief refine elements
        obj = refine(obj, refine_level);
        %> \brief draw mesh
        draw( obj, zvar );
        %> construction function
        function obj = NdgMesh1d( cell, Nv, vx, K, EToV, EToR )
            vy = zeros(size(vx)); % vy is all zeros
            vz = zeros(size(vx)); % vz is all zeros
            obj = obj@NdgMesh(cell, Nv, vx, vy, vz, K, EToV, EToR );
        end% func
        
        % assemble mesh connection
        function ConnectMeshUnion( obj, meshId, meshUnion )
            
        end
    end
    
end

