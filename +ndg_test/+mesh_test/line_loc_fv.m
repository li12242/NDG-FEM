classdef line_loc_fv < ndg_test.mesh_test.loc_fv
    %LINE_LOC_FV Summary of this class goes here
    %   Detailed explanation goes here
    
    methods
        function obj = line_loc_fv(mesh)
            obj = obj@ndg_test.mesh_test.loc_fv(mesh);
        end
    end
    
    %% Ë½ÓÐº¯Êý
    methods(Access=protected)
        function [ Nedge, Kloc, EToV ] = loc_connect(obj, cell) 
            Nedge = cell.N;
            Kloc = cell.N;
            EToV = zeros(2, Kloc);
            for k = 1:Kloc
                EToV(:, k) = [k, k+1]';
            end
        end% func
        
        function [ v1, v2, nx, ny, nz, ds ] = loc_edge_info(obj, mesh)
            K = mesh.K;
            v1 = zeros(obj.Nedge, K);
            v2 = zeros(obj.Nedge, K);
            nx = zeros(obj.Nedge, K);
            ny = zeros(obj.Nedge, K);
            nz = zeros(obj.Nedge, K);
            ds = ones(obj.Nedge, K);
            
            for e = 1:K
                sk = 1;
                Npoint = (e-1)*mesh.cell.Np;
                for k = 1:obj.Kloc
                    v1(sk, e) = obj.EToV(1,k)+Npoint;
                    v2(sk, e) = obj.EToV(2,k)+Npoint;
                    x1 = mesh.x( v1(sk, e) );
                    x2 = mesh.x( v2(sk, e) );
                    nx(sk, e) = sign(x2 - x1);
                end
            end
        end% func
        
        function [ vol ] = loc_fv_info(obj, mesh)
            Np = mesh.cell.Np;
            K = mesh.K;
            vol = zeros(Np, K);
            for e = 1:K
                for k = 1:obj.Kloc
                    v1 = obj.EToV(1, k);
                    v2 = obj.EToV(2, k);
                    d = abs( mesh.x(v1, e) - mesh.x(v2, e) );
                    vol(v1, e) = vol(v1, e) + d*0.5;
                    vol(v2, e) = vol(v2, e) + d*0.5;
                end
            end
        end% func
        
        function [P, R] = project_matrix(obj, mesh)
            cell = mesh.cell;
            P = zeros(cell.Np, cell.Np);
            
            [~, wq] = Polylib.zwgl(cell.N+1);
            Nq = numel(wq); wq = wq./sum(wq);
            Pw = zeros(Nq, cell.Np);
            len = zeros(cell.Np, 1);
            for k = 1:obj.Kloc
                v1 = obj.EToV(1, k);
                v2 = obj.EToV(2, k);
                r1 = cell.r(v1); 
                r2 = cell.r(v2);
                r0 = (r1 + r2)*0.5;
                
                rq1 = obj.project_vert2quad(cell.N, [r1, r0]');
                for i = 1:cell.Np
                    Pw(:,i) = cell.orthogonal_func(cell.N, i, rq1);
                end
                a = abs(r1 - r0);
                len(v1) = len(v1) + a;
                P(v1, :) = P(v1, :) + a*wq'*Pw;
                
                rq2 = obj.project_vert2quad(cell.N, [r0, r2]');
                for i = 1:cell.Np
                    Pw(:,i) = cell.orthogonal_func(cell.N, i, rq2);
                end
                a = abs(r2 - r0);
                len(v2) = len(v2) + a;
                P(v2, :) = P(v2, :) + a.*wq'*Pw;
            end
            P = bsxfun(@times, P, 1./len);
            P = P/cell.V;
            R = inv(P);
        end% func
        
        function quad_val = project_vert2quad(obj, N, vert_val)
            [rq, ~] = Polylib.zwgl(N+1);
            quad_val = ( (1-rq)*vert_val(1, :) ...
                + (1+rq)*vert_val(2, :) ).*0.5;
        end
    end
    
end

