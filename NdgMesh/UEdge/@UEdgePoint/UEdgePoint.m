classdef UEdgePoint < UEdge
    
    methods
        function obj = UEdgePoint(umesh, BCToV)
            obj = obj@UEdge(umesh, BCToV);
        end
    end
    
    methods(Hidden, Access=protected)
        
        function [bcell] = setStdEdgeCell(obj, N)
            bcell = StdPoint(0);
        end
        
        function [ bcType ] = assembleBC(obj, BCToV)
            tol = 1e-10;
            Nbc = size(BCToV, 2); % the column number
            bcType = ones(obj.Nedge, 1);
            for i = 1:Nbc
                for n = 1:obj.Nedge
                    if abs(BCToV(1, i) - obj.FToV(n)) < tol
                        bcType(n) = BCToV(2, i);
                    end
                end
            end
            bcType = EdgeBCType(bcType);
        end
        
        %> Calculate the outward normal vector and the determination of
        %> Jacobian matrix at quadrature points on each edge.
        function [ nxq, nyq, nzq, Jq ] = assembleQuadPointScale(obj, umesh)
            Nq = obj.bcell.Nq;
            nxq = ones(Nq, obj.Nedge);
            nyq = zeros(Nq, obj.Nedge);
            nzq = zeros(Nq, obj.Nedge);
            Jq = ones(Nq, obj.Nedge);
            
            xb1 = umesh.x( obj.FToN1 );
            % Define outward normals
            for n = 1:obj.Nedge
                k1 = obj.FToE(1, n);
                if abs(xb1(n) - min(umesh.x(k1)) ) < 1e-10
                    nxq(n) = -1;
                end
            end
        end
        
        function [ Nedge, FToE, FToF, FToV ] = assembleEdgeConnect( obj, umesh )
            [ EToFG ] = getUniqueFaceIndex(obj, umesh);
            [ EToE, EToF ] = assembleCellConnect(obj, umesh, EToFG);
            
            [~, id, ~] = unique(EToFG); % 寻找不同面标记号
            
            Nedge = numel(id);
            FToE = zeros(2, Nedge);
            FToF = zeros(2, Nedge);
            
            cell = umesh.cell;
            FToE(1, :) = fix( (id-1)./cell.Nface )+1; % local element index
            FToF(1, :) = rem(id-1, cell.Nface)+1; % 不同标记所在行号
            FToE(2, :) = EToE(id); % adjacent element index
            FToF(2, :) = EToF(id); % adjacent face index
            FToV = umesh.EToV(id);
        end% func
        
        
        function [EToE, EToF] = assembleCellConnect(obj, umesh, EToFG)
            cell = umesh.cell;
            EToE = ones(cell.Nface, 1)*(1:umesh.K);
            EToF = (1:cell.Nface)'*ones(1,umesh.K);
            
            for n = 1:(cell.Nface*umesh.K)
                m = find( abs(EToFG - EToFG(n))<1e-10 );
                t = m( m ~= n );
                if( ~isempty(t) )
                    EToE(n) = fix( (t-1)./cell.Nface )+1;
                    EToF(n) = rem(t-1, cell.Nface)+1;
                end
            end
        end
        
        function [ EToFG ] = getUniqueFaceIndex(obj, umesh)
            cell = umesh.cell;
            EToFG = zeros(cell.Nface, umesh.K);
            for f = 1:cell.Nface
                EToFG(f, :) = umesh.EToV(cell.FToV(1,f), :);
            end
        end% func
    end
    
end

