%> @brief Edge information on unstructured meshes.
%
%> And here we can put some more detailed informations about the class.
% ======================================================================
%> This class is part of the NDG-FEM software.
%> @author li12242, Tianjin University, li12242@tju.edu.cn
% ======================================================================
classdef UEdge < handle
    
    properties(Hidden=true, SetAccess=protected)
    end
    
    properties(SetAccess=protected)
        %> pointer to edge cell object
        bcell
        %> number of edges
        Nedge
        %> vertex index on each edge
        FToV
        %> adjacent cell index [kM; kP]
        FToE
        %> local face index of adjacent cells [fM; fP]
        FToF
        %> local node index of local cells
        FToN1
        %> local node index of adjacent cells
        FToN2
        %> boundary condition types of edges
        bcType @int8
        %> determination of Jacobian matrixs
        Jq
        %> x component of the outward normal vector
        nxq
        %> y component of the outward normal vector
        nyq
        %> z component of the outward normal vector
        nzq
    end
    
    methods(Abstract, Hidden=true, Access=protected)
        [ bcell ] = setStdEdgeCell(obj, N);
        [ bcType ] = assembleBC(obj, BCToV)
        [ Nedge, FToE, FToF, FToV ] = assembleEdgeConnect( obj, umesh );
        [ nxq, nyq, nzq, Jq ] = assembleQuadPointScale(obj, umesh);
    end
    
    methods
        function obj = UEdge(umesh, BCToV)
            
            [ obj.bcell ] = setStdEdgeCell( obj, umesh.cell.N );
            [ obj.Nedge, obj.FToE, obj.FToF, obj.FToV ] ...
                = assembleEdgeConnect( obj, umesh );
            
            [ obj.FToN1, obj.FToN2 ] = assembleNodeConnect( obj, umesh );
            [ obj.bcType ] = assembleBC(obj, BCToV);
            [ obj.nxq, obj.nyq, obj.nzq, obj.Jq ] ...
                = assembleQuadPointScale(obj, umesh);
            
        end
    end
    
    methods(Hidden, Access=protected)
        
            
%         function [ Nedge, FToV, FToE, FToF ] = assembleFaceConnect(obj, ecell, umesh)
%             
%             [ EToE ] = assembleCellConnect(obj, umesh);
%             [ EToF ] = assembleLocalFaceConnect(obj, umesh, EToE);
%             [Nedge, FToV, FToE, FToF] = assembleCellEdgeConnect(obj, EToV, EToE, EToF);
%         
%         end
%         
%         function [ EToE ] = assembleCellConnect(obj, umesh)
%             EToE = zeros(umesh.cell.Np, umesh.K);
%             for k = 1:umesh.K
%                 for f = 1:umesh.cell.Nface
%                     EToE(f, k) = getAdjacentCellIndex(...
%                         umesh.EToV, umesh.vx, umesh.vy, umesh.vz);
%                 end
%             end
%         end
        
%         function [ EToE, EToF ] = assembleCellConnect(obj, EToV)
%         end% func
%         
%         %> assemble the 
%         function [Nedge, FToV, FToE, FToF] = assembleCellEdgeConnect(obj, EToV, EToE, EToF)
%             ind = obj.mesh.Eind;
%             [~, id, ~] = unique(ind); % 寻找不同面标记号
%             
%             % 赋值
%             Nedge = numel(id);
%             FToE = zeros(2, Nedge);
%             FToF = zeros(2, Nedge);
%             FToV = zeros(2, Nedge);
%             
%             FToE(1, :) = fix( (id-1)./obj.cell.Nface )+1; % local element index
%             FToF(1, :) = rem(id-1, obj.cell.Nface)+1; % 不同标记所在行号
%             FToE(2, :) = EToE(id); % adjacent element index
%             FToF(2, :) = EToF(id); % adjacent face index
%             
%             FToV(1,:) = EToV(id);
%             FToV(2,:) = EToV(id);
%         end
        
%         function [ FToV ] = assembleNodeConnect(obj, Nedge, FToE, FToF)
%             
%         end
        function [ FToN1, FToN2 ] = assembleNodeConnect( obj, umesh )
        %function obj = node_connect(obj, mesh, Nedge, kM, kP, fM, fP)
            Nq = obj.bcell.Nq;
            FToN1 = zeros( Nq, obj.Nedge );
            FToN2 = zeros( Nq, obj.Nedge );
            
            Np = umesh.cell.Np;
            %sk = 1;
            for f = 1:obj.Nedge
                k1 = obj.FToE(1, f); 
                k2 = obj.FToE(2, f);
                f1 = obj.FToF(1, f); 
                f2 = obj.FToF(2, f);
%                 Nfp = obj.cell.Nfp(f1);
%                 list = 1:Nfp;
                
                nodeIndP = (k2-1)*Np + umesh.cell.Fmask(:, f2);
                for n = 1:Nq
                    FToN1(n, f) = (k1-1)*Np + umesh.cell.Fmask(n, f1);
                    %F(sk) = faceIndexStart(f1)+n;
                    xpM = umesh.x( FToN1(n, f) );
                    ypM = umesh.y( FToN1(n, f) );
                    zpM = umesh.z( FToN1(n, f) );
                    
                    xP = umesh.x( nodeIndP );
                    yP = umesh.y( nodeIndP );
                    zP = umesh.z( nodeIndP );
                    d12 = (xpM - xP).^2 + (ypM - yP).^2 + (zpM - zP).^2;
                    m = (d12 < 3e-16);
                    FToN2(n, f) = (k2-1)*Np + umesh.cell.Fmask(m, f2);
                end
            end% for
        end
        
        
    end
end

