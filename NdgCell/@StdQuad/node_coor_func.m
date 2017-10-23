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
function [ Np,r,s,t ] = node_coor_func( obj, order )
np = order+1;
[x,~] = zwglj(np);

r = x * ones(1, np);
s = ones(np, 1) * x';

r = r(:); 
s = s(:);
t = zeros(size(r));
Np = np*np;
end

