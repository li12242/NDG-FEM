function [ t, h ] = read_measured_data( filename )
%READ_MEASURED_DATA Summary of this function goes here
%   Detailed explanation goes here

fp = fopen(filename);
for i = 1:3
    str = fgetl(fp);
end
Np = fscanf(fp, '%d', 1);
data = fscanf(fp, '%f,%f\n', [2, Np]);
fclose(fp);

t = data(1, :);
h = data(2, :);
end

