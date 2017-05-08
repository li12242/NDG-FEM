classdef phys
    %PHYS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(Abstract)
        cell    % 标准单元对象
        mesh    % 网格对象
        Nfield  % 变量个数
        f_Q     % 变量
        f_extQ  % 外部值
    end
    
    methods(Abstract)
        disp(obj, field)
        vol_opt(obj, node_flux_func)
        suf_opt(obj, node_flux_func, nume_flux_func, obc_func)
        
    end
    
    methods
        function obj = phys(mesh, Nfield)
        end% func
        
        function obj = init_from_file(obj, filename)
            fp = fopen(filename);
            Num = fscanf(fp, '%d', 1);
            Nfld = fscanf(fp, '%d', 1); % read number of physical fields
            if ( ( (Num~=obj.mesh.K) && (Num~=obj.mesh.Nv) ) )
                error(['The number of values in file: ', ...
                    num2str(Num), ...
                    ' is neither element number: ', num2str(obj.mesh.K), ...
                    ' nor vertex number: ', num2str(obj.mesh.Nv)]);
            elseif (Nfld~=obj.Nfield)
                error(['The number of physical field in file: ', ...
                    num2str(Nfld), ...
                    ' is different from this phys object: ', ...
                    num2str(obj.Nfield)]);
            end
            fmtStr = ['%d ', repmat('%g ', 1, Nfld)];
            data = fscanf(fp, fmtStr, [Nfld+1, Num]);
            switch Num
                case obj.mesh.K
                    fprintf('\nInit with elemental averaged values.\n\n')
                    
                case obj.mesh.Nv
                    fprintf('\nInit with vertex values.\n\n')
            end
            fclose(fp);
        end
    end
    
end