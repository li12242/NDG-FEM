function ConvergenceRate
%% parameters
T        = 2.4;
ele      = [20, 40, 60, 80];
deg      = 1;
meshtype = 'tri';
filename = cell(numel(ele), 1);
for i = 1:numel(ele)
    filename{i} = ['Convection2D_', meshtype, '_', num2str(deg),'_',num2str(ele(i)), '.nc'];
end% for

%% postprocess class
PostproConv2d = Utilities.PostProcess.Postprocess(filename, meshtype, deg);

PrintTable    = PostproConv2d.GetConvTable('var', T, @ExactVar, 'M', ele')
end