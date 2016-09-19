%% Postprocess
% Class for post process
% 
% Meshods:
%   Postprocess     - construction function
%   GetVarData      - get the result of variable at spicific time
%   GetDofs         - get the number of unknows for 1D or 2D variables
%   NormErr         - calculate the norm error of variable
%   ConvRate        - calculate the convergence rate
%   Interp1d        - interpolation for 1-dimensional data
%   Interp2D        - interpolation for 2-dimensional data
%   SectProfile2D   - draw section profile for 2 dimension variables
%   Snapshot2D      - draw snapshot of 2 dimension variables
%   GetConvTable    - create table contains norm error and convergence rate
% 
% Warnning: 
%   The NetCDF files must contain 'x' or 'x' and 'y' as coordinate variables.
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
                obj.NcFile  = [obj.NcFile; ...
                    Utilities.PostProcess.ResultFile(filename{i})];
            end% for
            obj.StdCell = Utilities.PostProcess.StdCell(eleType, order);
        end% func
        
        %% Get convergence rate tables
        function PrintTable = GetConvTable(obj, varname, stime, exactFunH, varargin)
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
        
        %% Draw 1 dimensional snapshot of variables
        function p = Snapshot1D(obj, varname, stime, fileID, varargin)
            x   = obj.NcFile(fileID).GetVarData('x');
            [np, ne] = size(x);
            var      = obj.GetVarData(varname, stime, fileID);
            p = plot(x, var, varargin{:});
        end% func
        
        %% Draw 2 dimensional snapshot of variables
        function ph = Snapshot2D(obj, varname, stime, fileID)
            x   = obj.NcFile(fileID).GetVarData('x');
            y   = obj.NcFile(fileID).GetVarData('y');
            [np, ne] = size(x);
            var      = obj.GetVarData(varname, stime, fileID);
%             var(var<1e-2) = nan;
            vertex   = [x(:), y(:), var(:)];
            bclist   = obj.StdCell.bclist';
            EToV     = ones(ne, 1)*bclist;
            EToV     = EToV + (np*(0:ne-1))'*ones(size(bclist));
            ph = patch('Vertices', vertex, 'Faces', EToV,...
                'FaceColor', [0.8, 0.9, 1]);
        end% func
        
        function ph = SnapshotConst2D(obj, varname, fileID)
            x   = obj.NcFile(fileID).GetVarData('x');
            y   = obj.NcFile(fileID).GetVarData('y');
            [np, ne] = size(x);
            var      = obj.NcFile(fileID).GetVarData(varname);
            vertex   = [x(:), y(:), var(:)];
            bclist   = obj.StdCell.bclist';
            EToV     = ones(ne, 1)*bclist;
            EToV     = EToV + (np*(0:ne-1))'*ones(size(bclist));
            ph = patch('Vertices', vertex, 'Faces', EToV,...
                'FaceColor', [0.8, 0.9, 1]);
        end
        
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
        % Interp1d
        function yb = Interp1d(obj, numSol, xb, fileID)
            yb = zeros(size(xb));
            x      = obj.NcFile(fileID).GetVarData('x');
            for i = 1:numel(xb)
                yb(i) = disInterp1d(x, numSol, xb(i));
            end% for
        end% func
        
        % Interpolation for 2D variables, the variables shoule be arranged
        % the same with mesh coordinate
        function Data = Interp2D(obj, numSol, xc, yc, fileID)
            % check for input variable
            if isrow(xc)
                xc = xc';
            elseif isrow(yc)
                yc = yc';
            end% if
            % get position
            x      = obj.NcFile(fileID).GetVarData('x');
            y      = obj.NcFile(fileID).GetVarData('y');
            locVertList = obj.StdCell.verlist;
            % determine the cell inside
            flag   = Utilities.PostProcess.FindLocCell_Mex...
                (x,y,locVertList,xc,yc);
            
            % interpolation
            Data = zeros(size(xc));
            for i = 1:numel(xc)
                ind = find(flag == i);
                interp  = TriScatteredInterp...
                    (x(:,ind),y(:,ind),numSol(:,ind), 'linear');
                Data(i) = interp(xc,yc);
            end% for
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
            % get the output step, sk and skm1
            if (stime > time(end))
                error(['The input time %f is out of computation time range', ...
                    ' [%f, %f]\n'], stime, time(1), time(end));
            end% if
            sk = find(time>=stime, 1);
            if (sk==1)
                skm1 = 1;
            else 
                skm1 = sk - 1;
            end% if
            dt1 = abs(stime - time(skm1));
            dt2 = abs(stime - time(sk));
            w1  = dt2/(dt1+dt2);
            w2  = dt1/(dt1+dt2);
            % get the result
            sol1 = obj.NcFile(fileID).GetTimeVarData(varname, skm1);
            sol2 = obj.NcFile(fileID).GetTimeVarData(varname, sk);
            numSol = sol1*w1 + sol2*w2;
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
                    err = L2_2d(numSol, exSol);
%                     switch obj.StdCell.dim
%                         case 1
%                             err = L2_1d(numSol, exSol);
%                         case 2
%                             err = L2_2d(numSol, exSol);
%                     end% switch
                case 'Linf'
                    err = Linf_2d(numSol, exSol);
%                     switch obj.StdCell.dim
%                         case 1
%                             err = Linf_1d(numSol, exSol);
%                         case 2
%                             err = Linf_2d(numSol, exSol);
%                     end% switch
            end% switch
        end% func
    end% methods
end% classdef

%% disInterp
% Interpolation function for single node
function yb = disInterp1d(x, y, xb)
TOTAL = 1e-2;
% the boundary point is on vertex
ind   = abs(x-xb)<TOTAL;
if any(any( ind ))
    yb = mean(y(ind));
    return;
end% if
dx  = (x(1, :) - xb).*(x(end, :) - xb);
ind = find(dx < 0);
% interpolation with linear reconstruction
coef1 = (x(end, ind) - xb)/(x(end, ind) - x(1  , ind));
coef2 = (x(1  , ind) - xb)/(x(1  , ind) - x(end, ind));

yb    = y(1, ind)*coef1 + y(end, ind)*coef2;
end% func

%% Subroutine function
% functions for 
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