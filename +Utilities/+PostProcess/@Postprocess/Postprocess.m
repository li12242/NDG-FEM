%% Postprocess
% Class for post process
% 
% Meshods:
%   Postprocess     - construction function
% 
% 1.Norm errors:
%   GetDofs         - get the number of unknows for 1D or 2D variables
%   GetVarData      - get the result of variable at spicific time
%   NormErr         - calculate the norm error of variable
%   ConvRate        - calculate the convergence rate
%   GetConvTable    - create table contains norm error and convergence rate
% 
% 2.Interpolation:
%   Interp1D        - interpolation for 1-dimensional data
%   Interp2D        - interpolation for 2-dimensional data
% 
% 3.Plots:
%   Snapshot1D      - 
%   SectProfile2D   - draw section profile for 2 dimension variables
%   Snapshot2D      - draw snapshot of 2 dimension variables
%   SnapshotConst2D - draw snapshot of constant variables
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
        %   Postpro = Postprocess()
        function obj = Postprocess(filename, eleType, order)
            filenum     = numel(filename);
            obj.nfiles  = filenum;
            % create an array of ResultFile objects
            obj.NcFile  = []; 
            for i = 1:filenum
                obj.NcFile  = [obj.NcFile; ...
                    Utilities.PostProcess.ResultFile(filename{i})];
            end% for
            obj.StdCell = Utilities.PostProcess.StdCell(eleType, order);
        end% func
        
        %% 1. Norm errors function class
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
                case 'Linf'
                    err = Linf_2d(numSol, exSol);
            end% switch
        end% func
        
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
        
        PrintTable = GetConvTable(obj, varname, stime, exactFunH, varargin);
        
        %% 2. Interpolation function class
        function yb = Interp1D(obj, numSol, xb, fileID)
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
            eleId   = Utilities.PostProcess.FindLocCell_Mex...
                (x,y,locVertList,xc,yc);
            
            % interpolation
            Data = zeros(size(xc));
            for i = 1:numel(xc)
                ind = eleId(i);
                interp  = TriScatteredInterp...
                    (x(:,ind),y(:,ind),numSol(:,ind), 'linear');
                Data(i) = interp(xc(i),yc(i));
            end% for
        end% func
        
        %% 3. Plot function class
        
        function p = Snapshot1D(obj, varname, stime, fileID, varargin)
            x   = obj.NcFile(fileID).GetVarData('x');
            %[np, ne] = size(x);
            var      = obj.GetVarData(varname, stime, fileID);
            p = plot(x, var, varargin{:});
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
        
        function ph = SectProfile2D(obj, varname, stime, x, y, argname, fileID)
        % Usages:
        %   y = linspace(-4e3, 4e3, 100);
        %   x = zeros(size(y));
        %   PostproQuad.SectProfile2D('h', 0, x, y, 'y', 1)
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