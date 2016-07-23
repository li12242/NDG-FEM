%% Postprocess
% Class for post process
% 
% Warnning: 
% The NetCDF files must contain 'x' or 'x' and 'y' for coordinate variables.
% 
classdef Postprocess < handle
    properties
        nfiles   % number of files
        StdCell  % standard element object
        NcFile   % result files object
    end% properties
    
    methods
        %% Construction function
        % Input:
        %   filename - cell array contains warious resolution results
        %   eleType  - element type, 
        %   order    - order of degress
        % Usages:
        %
        function obj = Postprocess(filename, eleType, order)
            filenum     = numel(filename);
            obj.nfiles  = filenum;
            obj.NcFile  = [];
            for i = 1:filenum
                obj.NcFile  = [obj.NcFile; Utilities.PostProcess.ResultFile(filename{i})];
            end% for
            obj.StdCell = Utilities.PostProcess.StdCell(eleType, order);
        end% func
        
        %% Draw 2 dimensional snapshot of variables
        function Snapshot2D(obj, varname, stime, fileID)
            x   = obj.NcFile(fileID).GetVarData('x');
            y   = obj.NcFile(fileID).GetVarData('y');
            [np, ne] = size(x);
            var      = obj.GetVarData(varname, stime, fileID);
            vertex   = [x(:), y(:), var(:)];
            bclist   = obj.StdCell.bclist';
            EToV     = ones(ne, 1)*bclist;
            EToV     = EToV + (np*(0:ne-1))'*ones(size(bclist));
            patch('Vertices', vertex, 'Faces', EToV, 'FaceColor', [0.8, 0.9, 1])
        end% func
        
        %% Draw section profile for 2 dimensional variables
        % Usages:
        %
        %   y = linspace(-4e3, 4e3, 100);
        %   x = zeros(size(y));
        %   PostproQuad.SectProfile2D('h', 0, x, y, 'y', 1)
        %
        function ph = SectProfile2D(obj, varname, stime, x, y, argname, fileID)
            if nargin < 6
                error('Not enough input variables')
            end
            numSol = obj.GetVarData(varname, stime, fileID);
            data   = obj.Interp2D(numSol, x, y, fileID);
            switch argname
                case 'x'
                    ph = plot(x, data);
                case 'y'
                    ph = plot(y, data);
            end% switch
        end% func
        
        %% Interpolation
        % Interpolation for 2D variables, the variables shoule be arranged
        % the same with mesh coordinate
        function Data = Interp2D(obj, numSol, xc, yc, fileID)
            % check for variable dimensions
            x      = obj.NcFile(fileID).GetVarData('x');
            y      = obj.NcFile(fileID).GetVarData('y');
            Interp = scatteredInterpolant(x(:),y(:),numSol(:), 'linear');
            Data   = Interp(xc, yc);
        end% func
        
        %% Get the convergence rate
        function rate = ConvRate(obj, varname, stime, exFunH, errType)
            filenum = obj.nfiles;
            rate    = zeros(filenum, 1);
            err     = zeros(filenum, 1);
            % get error
            for i = 1:filenum
                err(i) = obj.NormErr( varname, stime, exFunH, errType, i);
            end% for
            % get dofs
            dofs    = obj.GetDofs;
            % get convergence rate
            for i = 2:filenum
                rate(i) = -log2( err(i-1)./err(i) )/log2( dofs(i-1)./dofs(i) );
            end
        end% func
        
        %% Get the number of unknows
        function dofs = GetDofs(obj)
            dofs  = zeros(obj.nfiles, 1);
            for i = 1:obj.nfiles
                switch obj.StdCell.dim
                    case 1 % for 1D problems, the DOFs is the points number
                        x = obj.NcFile(i).GetVarData('x');
                        dofs(i) = numel(x);
                    case 2 % for 2D problems, the DOFs is root of points number 
                        x = obj.NcFile(i).GetVarData('x');
                        dofs(i) = sqrt(numel(x));
                end% switch
            end% for
        end% func
        
        %% Get the result of variable
        function numSol = GetVarData(obj, varname, stime, fileID)
            time  = obj.NcFile(fileID).GetVarData('time');
            % get the output step
            [offsettime, ist] = min(abs(time - stime));
            if offsettime/stime > 0.1
                warning(['Attention: the time divergence is ', ...
                    num2str(offsettime)]);
            end% if
            % get the result
            numSol = obj.NcFile(fileID).GetTimeVarData(varname, ist);
        end% func
        
        %% Get the norm error
        function err = NormErr(obj, varname, stime, exFunH, errType, fileID)
            numSol = obj.GetVarData(varname, stime, fileID);
            % get the exact solution
            switch obj.StdCell.dim
                case 1
                    x = obj.NcFile(fileID).GetVarData('x');
                    exSol = exFunH(x, stime);
                case 2
                    x = obj.NcFile(fileID).GetVarData('x');
                    y = obj.NcFile(fileID).GetVarData('y');
                    exSol = exFunH(x, y, stime);
            end
            
            % get the Norm error
            switch errType
                case 'L2'
                    switch obj.StdCell.dim
                        case 1
                            err = L2_1d(numSol, exSol);
                        case 2
                            err = L2_2d(numSol, exSol);
                    end% switch
                case 'Linf'
                    switch obj.StdCell.dim
                        case 1
                            err = Linf_1d(numSol, exSol);
                        case 2
                            err = Linf_2d(numSol, exSol);
                    end% switch
            end% switch
        end% func
    end% methods
end% classdef

%% subroutine function
function err = L2_1d(y, ey)
err = sqrt( sum((y - ey).^2)./numel(y) );
end% func

function err = Linf_1d(y, ey)
err = max( abs(y - ey) );
end

function err = L2_2d(y, ey)
err = sqrt( sum2((y - ey).^2)./numel(y) );
end% func

function err = Linf_2d(y, ey)
err = max2( abs(y - ey) );
end

function sumM = sum2(M)
sumM = sum(sum(M));
end

function maxM = max2(M)
maxM = max(max(M));
end% func