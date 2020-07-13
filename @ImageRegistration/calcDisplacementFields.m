
function [fieldsX, fieldsY, error] = calcDisplacementFields (this, target_grads)
    fieldsX = {};
    fieldsY = {};
    error = [];
    
    for i = 1:size(target_grads,3)
        cuda_target_grad = gpuArray(target_grads(:,:,i));
        
        [intX, intY, intError] = this.calcIntegerShift(this.cPaddedTemplateGrad, cuda_target_grad, this.cPaddedPoints);
        
        valid_indx = find(intError <= (this.CCSize / 4));
        cuda_dispX = int32(intX(valid_indx));
        cuda_dispY = int32(intY(valid_indx));
        cuda_points = this.cPaddedPoints(valid_indx,:);
        
        [subX, subY] = this.calcSubpixelShift (this.cPaddedTemplateGrad, cuda_target_grad, cuda_points, cuda_dispX, cuda_dispY);
        
        translationX = double(cuda_dispX) + subX;
        translationY = double(cuda_dispY) + subY;
        
        translationX = gather(translationX);
        translationY = gather(translationY);
        points = double(this.Points(valid_indx,:));
        
        Fx = scatteredInterpolant(points(:,1), points(:,2), translationX, 'natural');
        Fy = scatteredInterpolant(points(:,1), points(:,2), translationY, 'natural');
        
        fieldsX{end + 1} = Fx;
        fieldsY{end + 1} = Fy;
        error(end + 1) = sqrt(mean(gather(intError) .^ 2));
        
        if (i == 1)
            this.cNextPaddedTemplate = cuda_target_grad;
        end
    end
    
    this.cPaddedTemplateGrad = this.cNextPaddedTemplate;
end