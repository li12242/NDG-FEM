%======================================================================
%> @brief Brief description of the function
%>
%> More detailed description.
%>
%> @param arg1 First argument
%> @param arg2 Second argument
%>
%> @retval out1 return value for the first output variable
%> @retval out2 return value for the second output variable
%======================================================================
%> This function is part of the NDGOM software. 
%> @author li12242, Tianjin University, li12242@tju.edu.cn
%======================================================================
function [ EToV ] = makeCounterclockwiseVertexOrder( EToV, vx, vy )
K = size(EToV, 2);
for k = 1:K
    vertId = EToV(:, k);
    vxk = vx( EToV(:, k) );
    vyk = vy( EToV(:, k) );

    vertOrder = convhull(vxk, vyk);
    EToV(:, k) = vertId( vertOrder(1:end-1) );
end

end

