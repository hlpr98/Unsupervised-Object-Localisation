addpath('./object_proposal/matlab');
addpath('./object_proposal/cmex');

%-------------- Global variables (needed for callbacks) --------------
global h imgId files imgDir configFile proposals;

%-------------- Input --------------
imgDir = 'test_images';
imgId = 1;
configFile = fullfile('./object_proposal/config', 'rp.mat');
%'config/rp_4segs.mat' to sample from 4 segmentations (slower but higher recall)
%'config/rp.mat' to sample from 1 segmentations (faster but lower recall)

%-------------- Find images in dir: --------------
files = dir(imgDir);
assert(numel(files) >= 3);
files = files(3 : end);
if(strcmp(files(1).name, '.svn'))
  files = files(2 : end);
end


%-------------- Processing: --------------
for i = 1:length(files)
    imgId = i;
    InteractiveCenterDemo(configFile);
end
