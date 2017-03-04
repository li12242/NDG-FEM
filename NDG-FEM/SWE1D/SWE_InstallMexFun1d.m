%% Install of mex function
oripath = pwd;
localpath = [oripath, '/SWE1D/Mex/'];
filename = dir(localpath);
cd(localpath)
for i = 1:numel(filename)
    if ~filename(i).isdir
        mex('-O', filename(i).name)
    end
end
cd(oripath)