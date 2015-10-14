%% varGroup
% 储存变量组
% 
classdef varGroup < handle
    properties(SetAccess = private)
        struc   % struct var to store variable
    end% properties
    
    methods
%% varGroup
% constructor function
        function obj = varGroup
        end% func
        
%% incert 
% 向 group 对象中添加变量
% INPUT:
%   varName     - 变量名称
%   varData     - 变量值
        function obj = incert(obj, varName, varData)
            obj.struc.(varName) = varData;
        end% func
%% getVal
% 获取 group 对象中储存的变量
% INPUT:
%   varName     - 变量名称
% OUTPUT:
%   data        - 变量值
        function data = getVal(obj, varName)
            data = obj.struc.(varName);
        end% func
%% print
% show the variable stored in the group
        function print(obj)
            fieldnames(obj.struc)
        end% func
    end% method
end% class