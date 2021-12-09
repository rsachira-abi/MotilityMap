
% Unit tests for ImageRegistration class

clear;
clc;

addpath('..');

%% Integer translation
disp('Running test: Integer translation');

[template, target] = Utils.generateShiftedPair(5, 7);

[meshX, meshY] = meshgrid(70:30:size(template,2) - 70, 70:30:size(template,1) - 70);
points = int32([meshX(:), meshY(:)]);

reg = ImageRegistration(template, points, 64);
[gradTemp, points] = reg.padArray(reg.grad(template), points);
[gradTarg, ~] = reg.padArray(reg.grad(target), points);

[intX, intY, intError] = reg.calcIntegerShift(gradTemp, gradTarg, points);

error = sqrt(mean((intX - 5).^2 + (intY - 7).^2))
try
    assert(error < 1e-2);
catch
    warning('Test Failed!');
    keyboard;
end

%% Subpixel translation test
disp('Running test: Subpixel translation');

[template, target] = Utils.generateShiftedPair(0.6, 0.7);

[meshX, meshY] = meshgrid(70:30:size(template,2) - 70, 70:30:size(template,1) - 70);
points = int32([meshX(:), meshY(:)]);

Disp = zeros(size(points,1), 1);

reg = ImageRegistration(template, points, 64);
[gradTemp, points] = reg.padArray(reg.grad(template), points);
[gradTarg, ~] = reg.padArray(reg.grad(target), points);

[subX, subY] = reg.calcSubpixelShift(gradTemp, gradTarg, points, gpuArray(Disp), gpuArray(Disp));

error = sqrt(mean((subX - 0.6).^2 + (subY - 0.7).^2))
try
    assert(error < 1e-2);
catch
    warning('Test Failed!');
    keyboard;
end

%% Translational shift
disp('Running test: Translation shift');

tr = 0.5:1.5:20;
LookAhead = 2;

Fields = cell(length(tr) + 1, 2);

[template, target] = Utils.generateShiftedPair(tr(1), 0);

[meshX, meshY] = meshgrid(70:30:size(template,2) - 70, 70:30:size(template,1) - 70);
points = int32([meshX(:), meshY(:)]);

reg = ImageRegistration(template, points, 64);

sz = reg.getTemplateSize();
Queue = complex(double.empty(sz(1), sz(2), 0));
qLength = 0;
TemplateIndx = 0;
for i = 1:length(tr)
    disp(['tr = ', num2str(i)]);
    
    [~, target] = Utils.generateShiftedPair(tr(i), 0);
    target_grad = reg.grad(target);
    [padded_target_grad, ~] = reg.padArray(target_grad, []);
    
    if (qLength < LookAhead)
        qLength = qLength + 1;
        Queue(:,:,qLength) = padded_target_grad;
        disp(['--- Added, length = ', num2str(qLength)]);
    end
    
    if (qLength == LookAhead)
        disp('--- Processing');
        [Fx, Fy] = reg.calcDisplacementFields(Queue);
        
        Fields{TemplateIndx + 1, 1} = Fx;
        Fields{TemplateIndx + 1, 2} = Fy;
        TemplateIndx = TemplateIndx + 1;
        
        Queue(:,:,1) = [];
        qLength = qLength - 1;
        disp(['--- Removed, length = ', num2str(qLength)]);
    end
end

while (qLength > 0)
    [Fx, Fy] = reg.calcDisplacementFields(Queue);
    
    Fields{TemplateIndx + 1, 1} = Fx;
    Fields{TemplateIndx + 1, 2} = Fy;
    TemplateIndx = TemplateIndx + 1;
    
    Queue(:,:,1) = [];
    qLength = qLength - 1;
    disp(['--- Removed, length = ', num2str(qLength)]);
    disp('--- Processing');
end

% Check the displacement values
points = double(points);
tr = [0, tr];
for i = 1:size(Fields,1)
    Fx = Fields{i,1};
    Fy = Fields{i,2};
    
    for j = 1:length(Fx)
        Fx_ = Fx{j};
        Fy_ = Fy{j};
        
        knownTranslation = tr(i + j) - tr(i);
        
        estimatedTranslation = Fx_(points(:,1), points(:,2));
        
        error = sqrt(mean((estimatedTranslation - knownTranslation).^2))
        try
            assert(error < 1e-2);
        catch
            warning('Test Failed!');
            keyboard;
        end
    end
end


















































