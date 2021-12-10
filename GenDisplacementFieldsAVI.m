%--------------------------------------------------------------------------
% GenDisplacementFieldsAVI:
%
% Generate the displacement fields from an AVI file between given frame
% numbers (inclusive).
%--------------------------------------------------------------------------

USE_CUDA = true;
IS_RAW = false;

project_dir = input('Project folder name: ');
[vFile, vPath] = uigetfile('*.avi', 'Select the video file');

%% Process video

StartFrameNum = input('Enter start frame number (default = 1): ');
EndFrameNum = input('Enter the end frame number (default = -1 for end): ');

vin = VideoReader(fullfile(vPath, vFile));

if (EndFrameNum == -1)
    EndFrameNum = vin.NumFrames;
end

NumFrames = EndFrameNum - StartFrameNum + 1;

FirstFrame = vin.read(StartFrameNum);
[FirstFrame, rect] = imcrop(FirstFrame);
FirstFrame = im2double(rgb2gray(FirstFrame));

%% Generate displacement fields

Width = size(FirstFrame, 2);
Height = size(FirstFrame, 1);

[MeshX, MeshY] = meshgrid(128:30:(Width - 128), 128:30:(Height - 128));
PointGrid = [MeshX(:), MeshY(:)];

if (USE_CUDA == true)
    reg = ImageRegistrationCUDA(FirstFrame, PointGrid, 64);
else
    reg = ImageRegistration(FirstFrame, PointGrid, 64);
end

DisplacementFields = cell(NumFrames, 4);

TemplateIndx = 0;
max_weight = 0;

for i = StartFrameNum:EndFrameNum
    disp(['FrameNum = ', i]);
    
    target = vin.read(i);
    target = im2double(rgb2gray(imcrop(target, rect)));
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
clear vin Width Height MeshX MeshY PointGrid reg TemplateIndx target padded_target_grad Fx Fy error weight

mkdir(project_dir);
save(fullfile(project_dir, 'displacement_fields.mat'), '-v7.3');
