classdef ugrid_map < handle
    %UGRID This class represent a map in ugrid
    properties
        set_source % Source set
        set_target % Target set
        source_dim  % dimension of source set
        target_dim  % dimension of target set
        size  % non-zero number of mappings
        imap  % mapping table in sparse matrix format
        name  % map name
    end
    
    methods
        function obj = ugrid_map(set_source, set_target, size, imap, name)
            %UGRID Construct an instance of this class
            obj.set_source = set_source;
            obj.set_target = set_target;
            obj.imap = imap;
            obj.size = size;
            obj.name = name;
        end
    end

end