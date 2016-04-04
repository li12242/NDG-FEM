function [x, eta] = ExactTranscriticalFlowData
% Output:
%   x - coordinate
%   eta - water level
chain = [398, 143; 418, 154; 432, 167; 448, 191;
    459, 214; 465, 229; 474, 255; 484, 292; 488, 309;
    493, 329; 497, 352; 502, 375; 505, 393; 509, 418; 515, 457];

leftPart = [127, 143; 398, 143];
rightPart = [516, 246; 966 246];

line = [leftPart; chain; rightPart];

[x, eta] = trans2xy(line(:, 1), line(:,2));
end


function [x, y] = trans2xy(pix_x, pix_y)

xmin = 0; xmax = 25;
ymin = 0; ymax = 0.5;

pix_xmin = 127; pix_xmax = 968;
pix_yup = 37; pix_ydown = 653;

x = (pix_x - pix_xmin)./(pix_xmax - pix_xmin).*(xmax - xmin);
y = 0.5 - (pix_y - pix_yup)./(pix_ydown - pix_yup).*(ymax - ymin);
end