classdef std_fv_cell
    %STD_FV_CELL Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(SetAccess=protected)
        cell % std_cell object
        NFV % # of finite volume
        NCV % # of control volume, equal to Np
        Nedge % # of edges inside, equal to NFV*3
        EToV % # of vertices in each finite volume
        
        v1, v2 % adjacent nodes of each edge
        vol % volume of each CV
        P, R % project and reconstruct matrix from lagrange basis 
             % coefficients to coutrol volume values, satisfying P*R = I
        LIFT
    end
    
    methods(Abstract, Access=protected, Hidden=true)
        [ NFV, EToV ] = loc_fv_info( obj )
        [ Nedge, v1, v2, ds ] = loc_edge_info( obj );
        [ vol ] = loc_cv_info(obj)
        [ P, R ] = loc_proj_recon_mat(obj)
    end
    
    methods
        function obj = std_fv_cell(cell)
            obj.cell = cell;
            obj.NCV = cell.Np;
            [ obj.NFV, obj.EToV ] = obj.loc_fv_info();
            [ obj.vol ] = loc_cv_info(obj);
            [ obj.Nedge, obj.v1, obj.v2 ] = obj.loc_edge_info();
            [ obj.P, obj.R ] = obj.loc_proj_recon_mat();
            [ obj.LIFT ] = obj.lift_mat();
        end
        
        function [ fv_val ] = project_node2fv(obj, node_val)
            fv_val = obj.P * node_val;
        end
        
        function [ node_val ] = reconst_fv2node(obj, fv_val)
            node_val = obj.R * fv_val;
        end
            
    end
    
    methods
        function [ LIFT ] = lift_mat(obj)
            scell = obj.cell;
            LIFT = zeros(scell.Np, scell.Nfptotal);
            for n = 1:scell.Nfptotal
                row = scell.Fmask(n);
                LIFT(row, n) = 1.0;
            end
        end
    end
    
end

