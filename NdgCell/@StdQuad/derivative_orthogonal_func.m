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
function [ fdr, fds, fdt ] = derivative_orthogonal_func(obj, N, ind, r, s, t)

% transform the index to two indexes.
i = mod( ind - 1, N + 1 );
j = floor( (ind - 1) / (N + 1) );
% calculate the derivative basis function values.
fdr = GradJacobiP(r(:), 0, 0, i).*JacobiP(s(:), 0, 0, j);
fds = JacobiP(r(:), 0, 0, i).*GradJacobiP(s(:), 0, 0, j);
fdt = zeros(size(fdr));
end

