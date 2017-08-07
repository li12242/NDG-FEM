classdef loc_fv
    %LOC_FV Local finite volume information
    %   loc_fv gathers the FV information in each elements.
    %
    
    %%  Ù–‘
    properties(SetAccess=protected)
        mesh
        Nedge % number of common edges in each element
        Kloc % number of sub-cells in each refined element
        EToV % node connection in sub-cell
        
        %% infomations for edges in each elements
        v1, v2 % adjacent nodes of each edge
        nx, ny, nz % normal vector of all edges pointing from v1 to v2
        ds % length/area of each edge inside elements
        P, R % project and reconstruct matrix from lagrange basis 
             % coefficients to finite volume values, satisfying P*R = I
        %% information for FV in each elements
        vol % length/area/volume of each FV in elements
        %% surface edge length/area
        Js
        LIFT
    end
    
    %% private methods
    methods(Abstract, Access=protected)
        [ Nedge, Kloc, EToV ] = loc_connect(obj, cell) % refine information
        [ v1, v2, nx, ny, nz, ds ] = loc_edge_info(obj, mesh)
        [ vol ] = loc_fv_info(obj, mesh)
        [ Js ] = loc_fv_surface(obj, mesh)
        [ P, R ] = project_matrix(obj, mesh)
    end
    
    %% public methods
    methods
        function obj = loc_fv(mesh)
            cell = mesh.cell;
            obj.mesh = mesh;
            % node connection in sub-cell
            [ obj.Nedge, obj.Kloc, obj.EToV ] = loc_connect(obj, cell); 
            % edge information
            [ obj.v1, obj.v2, obj.nx, obj.ny, obj.nz, obj.ds ] ...
                = loc_edge_info(obj, mesh);
            % finite volume information
            [ obj.vol ] = loc_fv_info(obj, mesh);
            [ obj.Js ] = loc_fv_surface(obj, mesh);
            obj.LIFT = zeros(cell.Np, cell.Nfptotal);
            for n = 1:cell.Nfptotal
                row = cell.Fmask(n);
                obj.LIFT(row, n) = 1.0;
            end
            % project matrix
            [ obj.P, obj.R ] = project_matrix(obj, mesh);
        end
        
        function val = project_node2fv(obj, f_Q)
            % project the DG results to fv results
            val = obj.P*f_Q;
        end
        
        function f_Q = project_fv2node(obj, val)
            f_Q = obj.R*val;
        end
        
    end
    
end

