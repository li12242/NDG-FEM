%======================================================================
%> @brief set the open boundary types for solving physical problems.
%>
%> The boundary conditions is necessory for many problems,
%> some porblmes need some time-depended external inofrmation to drive
%> the movement in the computation domain, so the user has to set the
%> obcType as NdgBCType.File and give an input file (NetCDF format)
%> which contains the exteranl values.
%> The user can also use the analytical functions to generate the
%> exteranl information with the obcType setting as NdgBCType.Func.
%> For the problem which do not require any external infomation, the
%> procedures for acquiring the information can be ignore, and the user
%> have to set the obcType as NdgBCType.None.
%>
%> @param obj Phys class
%> @param obcType enumeration for NdgBCType
%> @param varargin optional input for requiring external information.
%>
%> if the obcType equals to NdgBCType.File, optional input is
%> ( obcIntervalType, obcIntervalValue, obcNetcdfFileName ).
%> if the obcType equals to NdgBCType.File, optional input is
%> ( obcIntervalType, obcIntervalValue ).
%> if the obcType equals to NdgBCType.None, no optional input is required.
%>
%>
%======================================================================
%> This function is part of the NDGOM software.
%> @author li12242, Tianjin University, li12242@tju.edu.cn
%======================================================================
function setBoundaryType(obj, obcType, varargin)
% read obcType
try
    obj.obcType = NdgBCType( obcType );
catch
    msgID = 'PhysUnion:setOBC';
    msgtext = ['The input obcType should be an',...
        ' enumeration from the NdgBCType class.'];
    ME = MException(msgID, msgtext);
    throw(ME);
end
% read obcIntervalType
switch obcType
    case NdgBCType.None
        return;
    case {NdgBCType.Function, NdgBCType.File}
        try
            obj.obcIntervalType = NdgIntervalType( varargin{1} );
        catch
            msgID = 'PhysUnion:setOBC';
            msgtext = ['The input obcIntervalType should be an',...
                ' enumeration from the NdgIntervalType class.'];
            ME = MException(msgID, msgtext);
            throw(ME);
        end
end
% read obcIntervalValue
switch obj.obcIntervalType
    case NdgIntervalType.Constant
        obj.obcTimeInterval = varargin{2};
    case NdgIntervalType.DeltaTime
        obj.obcTimeInterval = varargin{2};
    case NdgIntervalType.DeltaStep
        obj.obcStepInterval = varargin{2};
    otherwise
end
% read obcFileName
if (obcType == NdgBCType.File)
    obj.obcFileName = varargin{3};
    checkObcNetcdfFile( obj.obcFileName );
end
end% func

%> todo finish check open boundary condition
function checkObcNetcdfFile( filename )

end