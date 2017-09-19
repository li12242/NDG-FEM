%> @brief Unstructured triangular mesh class.
%
%> The mesh object can be created through three ways:
% 
%> @par (1)<b>file</b>: 
%> The mesh files including three files, which are 
%> named by the case name with the suffix '.node', '.ele' and '.edge'. 
%> In these files, the node coordinates, the vertices index in 
%> each element, and the boudnary conditions are assigned.
% 
%> @par (2)<b>variable</b>: 
%> The user can also construct the mesh object by
%> giving the mesh informations directly. In this way, the variables `Nv`,
%> and the coordinates of the vertices `vx`, `vy` and `vz`, the number of
%> elements `K`, and the infomation of each element `EToV`, `EToR`, `EToBS`
%> should be given as the input properly.
% 
%> @par (3)<b>uniform</b>: 
%> It also flexible to just given the informations about the structured 
%> mesh and get the mesh object of this uniform mesh. The user should 
%> assign the range of the rectangle domain with `xlim` and `ylim`, the
%> number of elements on each coordinates `Mx` and `My`, and the boundary
%> types for four boundaries (southern, northern, western and eastern)
%> 
% ======================================================================
%> This class is part of the NDG-FEM software. 
%> @author li12242, Tianjin University, li12242@tju.edu.cn
% ======================================================================
classdef UMeshTri < UMeshUnion2d
    
    methods
        function obj = UMeshTri(cell, varargin)
            % check input element
            if( ~checkInputCellIsTri(cell) )            
                msgID = 'UMeshTri:inputError';
                msgtext = 'The input cell should be StdTri object.';
                ME = MException(msgID, msgtext);
                throw(ME);
            end
            
            obj = obj@UMeshUnion2d(cell, varargin{:});
        end% func
        
        obj = refine(obj, multi_rate); % ¼ÓÃÜÍø¸ñ
    end% methods
end


function [ isTri ] = checkInputCellIsTri(cell)
isTri = isa(cell, 'StdTri');
end