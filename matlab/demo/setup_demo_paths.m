function paths = setup_demo_paths()
%SETUP_DEMO_PATHS Add required folders to the MATLAB path for demos
demoDir = fileparts(mfilename('fullpath'));
matlabDir = fileparts(demoDir);
addpath(matlabDir);
addpath(fullfile(matlabDir,'utils'));
addpath(fullfile(matlabDir,'gof'));
addpath(fullfile(matlabDir,'plot'));
paths.demoDir = demoDir;
paths.matlabDir = matlabDir;
paths.rootDir = fileparts(matlabDir);
end

