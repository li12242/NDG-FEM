function [ i, j ] = trans_ind( N, ind )
%> @brief Transfer the linear mode index on the triangle to (i,j).
%> @details The mode sequence on the standard triangle is
%>   i = 0, j = 0,1,2,...,N;
%>   i = 1, j = 0,1,2,...,N-1;
%>   ...
%>   i = N, j = 0;
%> @param[in] N  maximum order
%> @param[in] ind linear index of the orthogonal basis
% ======================================================================
%> This function is part of the NDG-FEM software.
% ======================================================================
sk = 1;
for i = 0:N
    for j = 0:(N-i)
        if (abs(sk-ind)<10e-4)
            return;
        end
        sk = sk+1;
    end
end% for
end% func
