% Copyright (C) 2015 Tianjin University.
% 
% Please see the AUTHORS file in the main source directory for a full list
% of copyright holders.
% 
% li12242
% Department of Civil Engineering 
% Tianjin University Tianjin
% 
% li12242@tju.edu.cn
% 
% This library is free software; you can redistribute it and/or
% modify it under the terms of the GNU Lesser General Public
% License as published by the Free Software Foundation,
% version 2.1 of the License.
% 
% This library is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% Lesser General Public License for more details.
% 
% You should have received a copy of the GNU Lesser General Public
% License along with this library; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
% USA
%==========================================================================
%
% DESCRIPTION
%   
%
% INPUT
%    xx = explains
%
% OUTPUT
%    xx = explains
%
% EXAMPLE USAGE
%    
%==========================================================================
classdef array < handle
% ARRAY: vector of length n
    properties
        val     % value
        num
    end
%     properties(Dependent = true, SetAccess = private)
%         
%     end
    
    methods
        function obj = array(N)
            if nargin > 0
                obj.val = zeros(N,1);
                obj.num = N;
            end
        end %function
        
        function set.val(obj, value)
            obj.val = value;
        end
        
        function num = get.num(obj)
            num = numel(obj.val);
        end
    end %methods
end %classdef

