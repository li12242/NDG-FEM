function PrintTable = GetConvTable(obj, varname, stime, exactFunH, varargin)
% 
% Usages:
%   meshtype = 'tri';
%   elenum   = [80, 100, 120, 160, 200];
%   PostproSWE2D = Utilities.PostProcess.Postprocess(filename, meshtype, 1);
%   head = 'M', elenum';
% 
%   PrintTable = PostproSWE2D.GetConvTable('h', T, ...
%                @ParabolicBowl2d_ExtDepth, 'M', elenum');
%   PrintTable
% 

if mod(numel(varargin),2)
    warning(['The number of your input is ', numel(varargin),...
        ',the variable names and datas do not match']);
    error('Please check your input')
end

% create table
PrintTable   = table;
% set filename in PrintTable from input varargin, such as
%
% varargin{1}       varargin{3}
%   -----             -----    
% varargin{2}       varargin{4}
%
std = 1;
for i = 1:numel(varargin)/2
    % use user input str as field name
    PrintTable.(varargin{std}) = varargin{std+1};
    std = std + 2;
end
% get DOFs
PrintTable.dofs = obj.GetDofs;
% get norm error
errL2        = zeros(obj.nfiles, 1);
errLinf      = zeros(obj.nfiles, 1);
for i =1:obj.nfiles
    errL2(i)   = obj.NormErr(varname, stime, exactFunH, 'L2', i);
    errLinf(i) = obj.NormErr(varname, stime, exactFunH, 'Linf', i);
end% for
% convergence rate for variable
a2   = obj.ConvRate(varname, stime, exactFunH, 'L2');
ainf = obj.ConvRate(varname, stime, exactFunH, 'Linf');

PrintTable.L2   = errL2;
PrintTable.a2   = a2;
PrintTable.Linf = errLinf;
PrintTable.ainf = ainf;
end% func