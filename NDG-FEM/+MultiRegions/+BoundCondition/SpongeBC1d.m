classdef SpongeBC1d < MultiRegions.BoundCondition.SpongeBC
    %SPONGEBC1D Summary of this class goes here
    %   Detailed explanation goes here
    
properties
    xb    % boundary position of sponge layer
    width % width of the sponge layer
end

methods
    function obj = SpongeBC1d(BCflag, fileName, xb, xe)
        obj = obj@MultiRegions.BoundCondition.SpongeBC...
            (BCflag, fileName);
        % range of sponge layer
        obj.xb = xb;
        obj.width = abs(xe-xb);
    end% func
    
    %% GetAbsorpCoeff
    % return the absorption coefficient of each nodes. For inner 
    % computation element, the coefficient is 0, but for the sponge layer, 
    % the coefficient is obtained from
    % $$ \sigma = \sigma_{max} \left( \frac{x - x_b}{D} \right)^2 $$
    % where $D$ is the width of sponge layer.
    function sigma = GetAbsorpCoeff(obj, mesh, dt)
        sigma_max = .9/dt;
        sigma = zeros(size(mesh.x));
        x_spe = mesh.x(:, obj.SpEToE);
        sigma_spe = sigma_max*((x_spe - obj.xb)/obj.width).^2;
        sigma(:, obj.SpEToE) = sigma_spe;
    end% func
    
    %% GetBC
    % return the boundary condition of each nodes. For inner computation
    % domain, the value is 0 for all nodes, but for sponge layer element,
    % the external solution $h_e$ is obtained from the nearst boundary
    % vertex.
    % As the function should return the external solution at a spicific
    % time, an interpolation on time dimension is required. The function
    % will get the two time step that closest to the input time and
    % interpolate these external solutions as the approximation external 
    % solution.
    function var = GetBC(obj, mesh, varname, time)
        % check input time
        if time > obj.time(end)
            error('The require BC time: %f is larger than that in BC files', time);
        end% if
        % find the time step ind1 & ind2 and get the interpolation
        % coefficients
        ind2 = find(obj.time>=time, 1);
        if ind2 ==1 
            ind1  = 1;
            coef1 = 1; coef2 = 0;
        else
            ind1 = ind2 - 1;
            time1 = obj.time(ind1); time2 = obj.time(ind2);
            coef1 = abs( (time - time2)/(time1 - time2) );
            coef2 = 1 - coef1;
        end% if
        % get the external solution on vertex
        varB1 = obj.BCfile.GetTimeVarData(varname, ind1);
        varB2 = obj.BCfile.GetTimeVarData(varname, ind2);
        var_spe  = varB1*coef1 + varB2*coef2;
        % get the external solution on all nodes
%         np   = mesh.Shape.nNode;
        var  = zeros(size(mesh.x));
%         var_spe = varB(obj.SpEToBV);
        var(:, obj.SpEToE) = var_spe;
    end% func
end
    
end

%% FindNearestNode
% return the indices of elements in vector b which is closest to each
% elements of vector a
% 
% Input:
%   a - [ra x 1] vector
%   b - [rb x 1] vector
% Output:
%   bindex - [ra x 1] integer vector
% Usages:
%   a      = rand(3,1); 
%   b      = rand(4, 1);
%   [bind] = FindNearestNode(a,b);
% 
function [bindex] = FindNearestNode(a,b)
ra = numel(a);
rb = numel(b);

% transfer a and b into vector form
if (~iscolumn(a))
    a = a';
elseif(~iscolumn(b))
    b = b';
end

temp     = abs(ones(ra,1)*b' - a*ones(1, rb));
[~, ind] = min(temp,[], 2);
bindex   = ind';
end% func