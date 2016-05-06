function [r, s] = getCoor(nOrder)
% get the coordinate of standard quadrilateral
% coordinate 'r' comes first

% get points distribution on a edge
np = nOrder+1;
[x,~] = Polylib.zwglj(np);

% coord 'r' counts fist; 's' is second
r = x*ones(1, np);
s = ones(np, 1)*x';

r = r(:); s = s(:);
end% func