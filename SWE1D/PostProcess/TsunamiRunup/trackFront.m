function [time, frontPosition] = trackFront

filename = 'SWE1D.nc';
time = ncread(filename, 'time');
% timestep = numel(time)
x = ncread(filename, 'x');
startStep = 1; itime = startStep;
frontPosition = zeros(size(time));

% get result
h = ncread(filename, 'h', [1, itime],[inf, 1]);

% get transition element index
transIndex = TransitionCellIdentify(h);
frontPosition(1) = x(transIndex);
% bottom topography

for itime = 1:numel(time)
    try 
        h = ncread(filename, 'h', [1, itime],[inf, 1]);
        transIndex = TransitionCellIdentify(h);
        frontPosition(itime) = x(transIndex);
        fprintf('Processing: %f ...\n', itime/numel(time))
    catch 
        keyboard
    end% try
end% for

end% func

function transIndex = TransitionCellIdentify(h)
% identify the wet/dey interface element
% Input: 
%   h - water depth
%   bedElva - bottom elevation
% Output:
%   transIndex - bool flag for wet/dry transition element, size [1, Ne]

% identify transitation element
hPositive = 10^-3;
% define wet cells
iswet = (h > hPositive);

% when adjacent element possess different wet/dry status
% transIndex assignment is true
dsta = diff(iswet);
transIndex = find(dsta, 1);
% transitation element shoule be wet
transIndex = transIndex + 1;
end% func