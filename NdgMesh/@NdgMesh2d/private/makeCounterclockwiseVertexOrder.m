function [ EToV ] = makeCounterclockwiseVertexOrder( EToV, vx, vy )
K = size(EToV, 2);
for k = 1:K
    vertId = EToV(:, k);
    vxk = vx( EToV(:, k) );
    vyk = vy( EToV(:, k) );

    vertOrder = convhull(vxk, vyk);
    EToV(:, k) = vertId( vertOrder(1:end-1) );
end

end% func