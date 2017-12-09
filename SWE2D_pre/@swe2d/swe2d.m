classdef swe2d < ndg_lib.phys.phys2d
    %SWE2D 2d shallow water equations solver.
    %   The conservative from of the governing equations are 
    % 
    %   $$$$
    %   
    %   where $\mathbf{U} = (h, q_x, q_y)$ is the convervative variable, 
    %   $qx, qy$ is the flux rate along the x and y axis.
    
    properties(Constant)
        Nfield = 3 % 物理场个数，分别为 h，qx=hu，qy=hv
        gra = 9.80616 % 重力加速度
    end
    
    properties(Abstract, Constant)
        hmin
    end
    
    properties(SetAccess=protected)
        bot % bottom elevation
        bx, by % gradient of the topography terms
        ftime % final time
        dt % time interval
        wetflag % the wet element flag
        slopelimiter % slope limiter
    end
    
    %% private methods
    methods(Access=protected)
        function spe = char_len(obj, f_Q)
            % the character length of each nodes
            h = f_Q(:,:,1);
            q = sqrt( f_Q(:,:,2).^2 + f_Q(:,:,3).^2 );
            spe = (q./h) + sqrt(obj.gra*h);
            spe(:, ~obj.wetflag) = eps;
        end
        
        function wetdry_detector(obj, f_Q)
            % get the wet/dry status for each elements
            
            % The elements with the the water depth of all the nodes
            % are greater than the threshold are set as the wet elements
            obj.wetflag = all( f_Q(:,:,1) > obj.hmin ); 
            obj.mesh.EToR( ~obj.wetflag ) = ndg_lib.mesh_type.Dry;
            obj.mesh.EToR( obj.wetflag ) = ndg_lib.mesh_type.Normal;
        end
        
        function topo_grad_term(obj)
            % Pre-calculate the gradient term of the bottom topography.
            obj.bx = obj.mesh.rx.*(obj.mesh.cell.Dr*obj.bot) ...
               + obj.mesh.sx.*(obj.mesh.cell.Ds*obj.bot);
            obj.by = obj.mesh.ry.*(obj.mesh.cell.Dr*obj.bot) ...
               + obj.mesh.sy.*(obj.mesh.cell.Ds*obj.bot);
        end
        
        [ rhs ] = rhs_term(obj, f_Q ) % Calculate the R.H.S. term
        [ dflux ] = surf_term(obj, f_Q ) % Calculate the surface flux term
        [ E,G ] = flux_term(obj, f_Q ) % Calculate the volume flux term
        [ sb ] = topo_sour_term( obj, f_Q ); % the topography source term
        [ dflux ] = hll_surf_term(obj, f_Q ) % HLL numerical flux term
        [ dflux ] = lf_surf_term(obj, f_Q ) % LF numerical flux term
    end
    
    %% Public methods
    methods(Abstract)
        [ f_Q ] = init(obj)
    end
    
    methods
        [ obj ] = RK45(obj); % 采用 SSP RK45 格式计算
        [ obj ] = RK45_OBC(obj); % 采用 SSP RK45 格式计算，并考虑开边界条件
        [ obj ] = VB_RK45(obj); % 采用 vertex-based RK 格式计算
        [ obj ] = VB_RK45_OBC(obj); % 采用 VB-RK 格式计算，并考虑开边界条件
        
        function obj = swe2d(varargin)
            switch nargin
                case 1
                    mesh = varargin{1};
                case 3
                    mesh = read_from_file(varargin{:});
                case 7
                    check_input(varargin{:});
                    mesh = uniform_mesh(varargin{:});
            end
            obj = obj@ndg_lib.phys.phys2d(mesh);
            
            %;
        end
    end
    
end

function check_input(varargin)
% N = varargin{1};
xlim = varargin{2};
ylim = varargin{3};
if ( xlim(2) < xlim(1) ) || ( ylim(2) < ylim(1) )
    error('The input "xlim" and "ylim" should be increasing!');
end% if
% Mx = varargin{4};
% My = varargin{5};
cell_type = varargin{6};
if ( ~isa(cell_type, 'ndg_lib.std_cell_type') )
    error('The input "cell_type" is not "ndg_lib.std_cell_type" object!');
end
face_type = varargin{7};
if ( ~isa(face_type, 'ndg_lib.bc_type') )
    error('The input "cell_type" is not "ndg_lib.std_cell_type" object!');
end
end% func

function mesh = read_from_file(varargin)
N = varargin{1};
casename = varargin{2};
cell_type = varargin{3};
switch cell_type
    case ndg_lib.std_cell_type.Tri
        cell = ndg_lib.std_cell.tri(N);
        mesh = ndg_lib.mesh.tri_mesh(cell, 'file', casename);
    case ndg_lib.std_cell_type.Quad
        cell = ndg_lib.std_cell.quad(N);
        mesh = ndg_lib.mesh.quad_mesh(cell, 'file', casename);
end% switch
end% func

function mesh = uniform_mesh(varargin)
N = varargin{1};
xlim = varargin{2};
ylim = varargin{3};
Mx = varargin{4};
My = varargin{5};
cell_type = varargin{6};
face_type = varargin{7};

switch cell_type
    case ndg_lib.std_cell_type.Tri
        cell = ndg_lib.get_std_cell(N, cell_type);
        mesh = ndg_lib.mesh.tri_mesh...
            (cell, 'uniform', {xlim, ylim, Mx, My, face_type});
        
    case ndg_lib.std_cell_type.Quad
        cell = ndg_lib.get_std_cell(N, cell_type);
        mesh = ndg_lib.mesh.quad_mesh...
            (cell, 'uniform', {xlim, ylim, Mx, My, face_type});
end
end% func

