function GenHeaderFile(n)
% Generate header files for DGOM model
% 
% Input
%   n : order of degree (less than 10!)
% Output:
%   file: 'dataN0X.h' (X is the order n)

shape = StdRegions.Quad(n);
filename = ['dataN0', num2str(n),'.h'];

% open file to write
fig = fopen(filename, 'w');
fprintf(fig, '#ifndef DATA0%d\n', n);
fprintf(fig, '#define DATA0%d 1\n\n', n);
% coordinate
format = '%18.14e,';
format = repmat(format, 1, shape.nNode);
format(end) = [];
fprintf(fig, ['double p_r[%d] = {', format,'};\n'], ...
    shape.nNode, shape.r);
fprintf(fig, ['double p_s[%d] = {', format,'};\n\n'], ...
    shape.nNode, shape.s);

% Integral coefficient
[~,w] = Polylib.zwglj(n+1);
format = '%18.14e,';
format = repmat(format, 1, n+1);
format(end) = [];
fprintf(fig, ['double p_w[%d] = {', format,'};\n\n'], ...
    n+1, w);

% Volume integral coeff
we = sum(shape.M);
format = '%18.14e,';
format = repmat(format, 1, shape.nNode);
format(end) = [];
fprintf(fig, ['double p_we[%d] = {', format,'};\n'], ...
    shape.nNode, we);

% derivative matrix
format = '%18.14e,';
format = repmat(format, 1, shape.nNode);
format(end) = [];
fprintf(fig, ['double p_Dr[%d][%d] = {{', format,'},\n'], ...
    shape.nNode, shape.nNode, shape.Dr(1, :));
for i = 2:shape.nNode-1
    fprintf(fig, ['{', format,'},'], shape.Dr(i, :));
end% for
fprintf(fig, ['{', format,'}};\n\n'], shape.Dr(end, :));

fprintf(fig, ['double p_Ds[%d][%d] = {{', format,'},\n'], ...
    shape.nNode, shape.nNode, shape.Ds(1, :));
for i = 2:shape.nNode-1
    fprintf(fig, ['{', format,'},'], shape.Ds(i, :));
end% for
fprintf(fig, ['{', format,'}};\n\n'], shape.Ds(end, :));

% LIFT matrix
nfp = size(shape.LIFT, 2); % number of points on all faces
format = '%18.14e,';
format = repmat(format, 1, nfp);
format(end) = [];

fprintf(fig, ['double p_LIFT[%d][%d] = {{', format,'},\n'], ...
    shape.nNode, nfp, shape.LIFT(1, :));
for i = 2:shape.nNode-1
    fprintf(fig, ['{', format,'},'], shape.LIFT(i, :));
end% for
fprintf(fig, ['{', format,'}};\n\n'], shape.LIFT(end, :));


% Fmask
nfp = nfp/shape.nFace;
fmask = zeros(shape.nFace, nfp);
for i = 1:shape.nFace
    [~, temp] = shape.getNodeListAtFace(i);
    fmask(i, :) = temp(:);
end%for
fmask = fmask -1; % index start from 0

format = '%d,';
format = repmat(format, 1, nfp);
format(end) = [];
fprintf(fig, ['int p_Fmask[%d][%d] = {{', format,'},\n'], ...
    shape.nFace, nfp, fmask(1, :));
for i = 2:shape.nFace-1
    fprintf(fig, ['{', format,'},'], fmask(i, :));
end% for
fprintf(fig, ['{', format,'}};\n\n'], fmask(end, :));


% finish
fprintf(fig, '#endif\n');
fclose(fig);
end% func