classdef UEdgeLine < UEdge
    %UEDGELINE Summary of this class goes here
    %   Detailed explanation goes here
    
    methods
        function obj = UEdgeLine(umesh, BCToV)
            obj = obj@UEdge(umesh, BCToV);
        end
    end
    
    methods(Hidden, Access=protected)
        function [bcell] = setStdEdgeCell(obj, N)
            bcell = StdLine(N);
        end
        
        function [ bcType ] = assembleBC(obj, BCToV)
            tol = 1e-10;
            Nbc = size(BCToV, 2); % the column number
            bcType = ones(obj.Nedge, 1);
            for i = 1:Nbc
                for n = 1:obj.Nedge
                    if (abs(BCToV(1, i) - obj.FToV(1, n)) < tol ) && ...
                        ( abs(BCToV(2, i) - obj.FToV(2, n)) < tol )
                        bcType(n) = BCToV(3, i);
                    end
                end
            end
            bcType = EdgeBCType(bcType); % transform to BoundaryType
        end
        
        function [ nxq, nyq, nzq, Jq ] = assembleQuadPointScale(obj, umesh)
            Nq = obj.bcell.Nq;
            nxq = zeros(Nq, obj.Nedge);
            nyq = zeros(Nq, obj.Nedge);
            nzq = zeros(Nq, obj.Nedge);
            Jq = zeros(Nq, obj.Nedge);
            
            for n = 1:obj.Nedge
                k1 = obj.FToE(1, n);
                f1 = obj.FToF(1, n);
                vert = umesh.EToV(umesh.cell.FToV(:, f1), k1);
                vx = umesh.vx(vert);
                vy = umesh.vy(vert);
                
                nx =  ( vy(2) - vy(1) );
                ny = -( vx(2) - vx(1) );
                
                J = sqrt(nx.*nx+ny.*ny); 
                nxq(:, n) = nx./J; 
                nyq(:, n) = ny./J;
                Jq(:, n) = J.*0.5;
            end
        end
        
        function [ Nedge, FToE, FToF, FToV ] = assembleEdgeConnect( obj, umesh )
            [ EToFG ] = getUniqueFaceIndex(obj, umesh);
            [ EToE, EToF ] = assembleCellConnect(obj, umesh, EToFG);
            
            [~, id, ~] = unique(EToFG); % 寻找不同面标记号
            
            Nedge = numel(id);
            FToE = zeros(2, Nedge);
            FToF = zeros(2, Nedge);
            FToV = zeros(2, Nedge);
            
            cell = umesh.cell;
            FToE(1, :) = fix( (id-1)./cell.Nface )+1; % local element index
            FToF(1, :) = rem(id-1, cell.Nface)+1; % 不同标记所在行号
            FToE(2, :) = EToE(id); % adjacent element index
            FToF(2, :) = EToF(id); % adjacent face index
            
            v1 = sub2ind([umesh.cell.Nv, umesh.K], ...
                umesh.cell.FToV(1, FToF(1, :)), FToE(1, :));
            v2 = sub2ind([umesh.cell.Nv, umesh.K], ...
                umesh.cell.FToV(2, FToF(1, :)), FToE(1, :));
            FToV(1, :) = umesh.EToV(v1);
            FToV(2, :) = umesh.EToV(v2);
            
        end
        
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
        end% func
        
        function [ EToFG ] = getUniqueFaceIndex(obj, umesh)
            cell = umesh.cell;
            EToFG = zeros(cell.Nface, umesh.K);
            for f = 1:cell.Nface
                v1 = umesh.EToV(cell.FToV(1,f), :);
                v2 = umesh.EToV(cell.FToV(2,f), :);
                EToFG(f, :) = min(v1, v2)*umesh.Nv + max(v1, v2);
            end            
        end% func
    end
    
end

