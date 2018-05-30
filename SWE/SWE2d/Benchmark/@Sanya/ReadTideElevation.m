function ReadTideElevation( obj )
%READTIDEELEVATION 此处显示有关此函数的摘要
%   此处显示详细说明
    [row, col] = size(obj.OBVid);
    
    Tide_file =[pwd, '/SWE2d/@Sanya/tide/TideElevation0112.txt'];
    fp = fopen(Tide_file,'r');
    obj.Tide = fscanf(fp, '%f\n', [row, inf]);
    fclose(fp);
%     obj.Tide = (obj.Tide)';

end

