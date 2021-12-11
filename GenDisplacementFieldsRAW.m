%--------------------------------------------------------------------------
% GenDisplacementFieldsRAW:
%
% Generate displacement fields from RAW image files between given frame
% numbers (inclusive).
%--------------------------------------------------------------------------

USE_CUDA = true;
IS_RAW = true;

project_dir = input('Project folder name: ', 's');
vPath = uigetdir('', 'Folder containing raw frames');

%% Process video

% Get all *.raw files
all_files = dir(fullfile(vPath, '*.raw'));
all_files = struct2cell(all_files);
all_files = cellfun(@extractFrameNumCellFunc , all_files(1,:));
all_files = sort(all_files);

StartFrameName = input(['Enter start frame name (', num2str(all_files(1)), '): ']);
EndFrameName = input(['Enter the end frame name (', num2str(all_files(end)), '): ']);

StartFrameNum = find(all_files == StartFrameName);
EndFrameNum = find(all_files == EndFrameName);
FrameNameList = all_files(StartFrameNum:EndFrameNum);
NumFrames = length(FrameNameList);

FirstFrame = ImportRAW(USE_CUDA, fullfile(vPath, [num2str(FrameNameList(1)), '.raw']));
[FirstFrame, rect] = imcrop(FirstFrame);

%% Generate displacement fields

Width = size(FirstFrame, 2);
Height = size(FirstFrame, 1);

[MeshX, MeshY] = meshgrid(128:30:(Width - 128), 128:30:(Height - 128));
PointGrid = int32([MeshX(:), MeshY(:)]);

if (USE_CUDA == true)
    reg = ImageRegistrationCUDA(FirstFrame, PointGrid, 64);
else
    reg = ImageRegistration(FirstFrame, PointGrid, 64);
end

DisplacementFields = cell(NumFrames, 4);

TemplateIndx = 0;
max_weight = 0;

for i = 1:NumFrames
    FrameName = [num2str(FrameNameList(i)), '.raw'];
    disp(['Frame = ', FrameName]);
    
    target = ImportRAW(USE_CUDA, fullfile(vPath, FrameName));
    target = imcrop(target, rect);
    target_grad = reg.grad(target);
    [padded_target_grad, ~] = reg.padArray(target_grad, []);
    
    [Fx, Fy, error] = reg.calcDisplacementFields(padded_target_grad);
    weight = 1 ./ error;
    weight(weight > 1) = 1;
    max_weight = max(weight);
    
    TemplateIndx = TemplateIndx + 1;
    DisplacementFields{TemplateIndx, 1} = Fx;
    DisplacementFields{TemplateIndx, 2} = Fy;
    DisplacementFields{TemplateIndx, 3} = i;
    DisplacementFields{TemplateIndx, 4} = weight;
end

% clear unwanted variables
clear StartFrameName EndFrameName Width Height MeshX MeshY PointGrid reg TemplateIndx target padded_target_grad Fx Fy error weight

mkdir(project_dir);
save(fullfile(project_dir, 'displacement_fields.mat'), '-v7.3');

%% Helper Function

function I = ImportRAW (USE_CUDA, file_path)
    if (USE_CUDA == true)
        I = ImageRegistrationCUDA.ImportRaw(file_path);
    else
        I = ImageRegistration.ImportRaw(file_path);
    end
end

function y = extractFrameNumCellFunc (x)
% extractFrameNumCellFunc
% Designed to be run with cellfun(). Extract the numerical frame number
% from the raw file listing names obtained with dir().
%
% Input:
%           x = Name of a single raw file. e.g.: '100.raw'
% Output:
%           y = Frame number extracted from the raw file name. e.g.: 100

    str_parts = split(x, '.');
    y = str2double(str_parts{1});
end



