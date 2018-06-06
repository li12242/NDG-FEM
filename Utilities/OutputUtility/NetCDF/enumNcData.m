classdef enumNcData < int8
    
    enumeration
        %> NAT = 'Not A Type' (c.f. NaN)
        NC_NAT  (0)
        %> signed 1 byte integer
        NC_BYTE (1)
        %> ISO/ASCII character
        NC_CHAR (2)
        %> signed 2 byte integer
        NC_SHORT (3)
        %> signed 4 byte integer
        NC_INT (4)
        %> single precision floating point number
        NC_FLOAT (5)
        %> double precision floating point number
        NC_DOUBLE (6)
    end
    
end

