function [EToE, EToF] = ele_connect(obj, EToV)
%ELE_CONNECT Summary of this function goes here
%   Detailed explanation goes here

% count the unique edges for triangle
ind = zeros(obj.cell.Nface, obj.K);
for f = 1:obj.cell.Nface
    % find two vertex index of each edge
    v1 = EToV(obj.cell.FToV(1,f), :);
    % calculate the indicator for each edge
    ind(f, :) = v1;
end

EToE = ones(obj.cell.Nface, 1)*(1:obj.K); % initialize for EToE
EToF = (1:obj.cell.Nface)'*ones(1,obj.K); % initialize for EToF

for n = 1:(obj.cell.Nface*obj.K)
    m = find( abs(ind - ind(n))<1e-10 );
    t = m( m ~= n );
    if( ~isempty(t) )
        EToE(n) = fix( (t-1)./obj.cell.Nface )+1;
        EToF(n) = rem(t-1, obj.cell.Nface)+1;
    end
end

end

