%======================================================================
%> @brief Rotate the node counter-clockwise with theta radian
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
function [ x, y ] = rotateNodeCoordinate2d( x, y, xc, yc, theta )

dx = x - xc;
dy = y - yc;

dx1 = cos(theta) * dx - sin(theta) * dy;
dy1 = sin(theta) * dx + cos(theta) * dy;

x = dx1 + xc;
y = dy1 + yc;

end% func

