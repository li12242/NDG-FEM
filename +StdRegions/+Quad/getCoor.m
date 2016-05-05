function [r, s] = getCoor(nOrder)
% get the coordinate of standard quadrilateral

% get points distribution on edges
np = nOrder+1;
[x,~] = Polylib.zwglj(np);

% coord 'r' fist; coord 's' second
r = x*ones(1, np);
s = ones(np, 1)*x';

r = r(:); s = s(:);
end% func