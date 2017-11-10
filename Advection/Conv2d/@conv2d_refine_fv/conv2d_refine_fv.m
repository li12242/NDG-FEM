classdef conv2d_refine_fv < conv2d_advection
    %CONV2D_REFINE_FV 2D convection problem with FVM for refined 
    %control volume
    %   Detailed explanation goes here
    
    properties(Constant)
        detect_threshold = 10;
    end
    
    properties
        loc_fv
        edge
    end
    
    methods
        [ obj ] = RK45( obj );
        
        function obj = conv2d_refine_fv(varargin)
            obj = obj@conv2d_advection(varargin{:});
            
            obj.edge = ndg_lib.mesh.edge(obj.mesh);
            stdcell_fv = ndg_test.cell_test.std_fv_tri(obj.mesh.cell);
            obj.loc_fv = ndg_test.mesh_test.mesh_loc_tri(stdcell_fv, obj.mesh);
        end
        
        function detect_refine_cell(obj)
%             res =  obj.mesh.rx.*(obj.mesh.cell.Dr*obj.f_Q) ...
%                 + obj.mesh.sx.*(obj.mesh.cell.Ds*obj.f_Q) ...
%                 + obj.mesh.ry.*(obj.mesh.cell.Dr*obj.f_Q) ...
%                 + obj.mesh.sy.*(obj.mesh.cell.Ds*obj.f_Q);
%             
%             refine_cell = any( abs(res) > obj.detect_threshold);
%             obj.mesh.EToR = ndg_lib.mesh_type.Normal;
%             obj.mesh.EToR( refine_cell ) = ndg_lib.mesh_type.Refine;
            
            obj.mesh.EToR(:) = ndg_lib.mesh_type.Refine;
        end
        
        function f_Q = refine2fv( obj, f_Q )
            v_Q = obj.loc_fv.project_node2fv( f_Q );
            ind = obj.mesh.EToR == ndg_lib.mesh_type.Refine;
            f_Q(:, ind) = v_Q(:, ind);
        end
        
        function f_Q = combine2node( obj, f_Q )
            v_Q = obj.loc_fv.project_fv2node(f_Q);
            ind = obj.mesh.EToR == ndg_lib.mesh_type.Refine;
            f_Q(:, ind) = v_Q(:, ind);
        end
    end
    
    methods(Access=protected) % private
        [ rhsQ ] = rhs_term( obj, f_Q ) % get the r.h.s term
    end
    
end
