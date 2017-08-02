classdef conv2d_refine_fv < conv2d_advection
    %CONV2D_REFINE_FV 2D convection problem with FVM for refined 
    %control volume
    %   Detailed explanation goes here
    
    properties(Constant)
        detect_threshold = 10;
    end
    
    properties
        v_Q % results in sub-cell control volume
    end
    
    methods
        [ obj ] = RK45( obj );
        
        function obj = conv2d_refine_fv(varargin)
            obj = obj@conv2d_advection(varargin{:});
        end
        
        function detect_refine_cell(obj)
            res =  obj.mesh.rx.*(obj.mesh.cell.Dr*obj.f_Q) ...
                + obj.mesh.sx.*(obj.mesh.cell.Ds*obj.f_Q) ...
                + obj.mesh.ry.*(obj.mesh.cell.Dr*obj.f_Q) ...
                + obj.mesh.sy.*(obj.mesh.cell.Ds*obj.f_Q);
            
            refine_cell = any( abs(res) > obj.detect_threshold);
            obj.mesh.EToR( refine_cell ) = ndg_lib.mesh_type.Refine;
            obj.mesh.EToR( ~refine_cell ) = ndg_lib.mesh_type.Normal;
        end
    end
    
    methods(Access=protected) % private
        [ E, G ] = flux_term( obj, f_Q ) % get the flux terms
        [ dflux ] = surf_term( obj, f_Q ) % get flux deviation
        [ rhs ] = rhs_term(obj, f_Q ) % get the r.h.s term
    end
    
end
