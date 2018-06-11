%> @brief emumeration for mesh region types.
%
%> And here we can put some more detailed informations about the class.
% ======================================================================
%> This class is part of the NDGOM software. 
%> @author li12242, Tianjin University, li12242@tju.edu.cn
% ======================================================================
classdef enumRegion < int8
    enumeration
        Normal      (0) % normal element
        Coarse      (1) % coarse element, exclude from calculation
        Refine      (2) % refined element, participate in calculation
    end
end