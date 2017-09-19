classdef mesh_loc_fv
    %MESH_LOC_FV Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        subcell
        mesh
        nx, ny, nz
        ds
        vol
        Js
    end
    
    methods(Access=protected, Hidden=true)
        [ nx, ny, nz, ds ] = inner_edge_info(obj)
        [ vol ] = loc_cv_vol(obj)
        [ Js ] = cv_surf_info(obj)
    end
    
    methods
        function obj = mesh_loc_fv(subcell, mesh)
            obj.subcell = subcell;
            obj.mesh = mesh;
            [ obj.nx, obj.ny, obj.nz, obj.ds ] = obj.inner_edge_info();
            [ obj.vol ] = obj.loc_cv_vol();
            [ obj.Js ] = obj.cv_surf_info();
        end% func
        
        function [fv_val] = project_node2fv(obj, node_val)
            fv_val = obj.subcell.P * node_val;
        end% func
        
        function [node_val] = project_fv2node(obj, fv_val)
            node_val = obj.subcell.R * fv_val;
        end% func
    end
    
end

