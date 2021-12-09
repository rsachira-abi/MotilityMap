function [subX, subY] = calcSubpixelShift (this, PaddedTemplateGrad, PaddedTargetGrad, Points, DispX, DispY)
    
    roi_1 = complex(zeros(this.CCSize, this.CCSize, size(Points,1)));
    roi_2 = complex(zeros(this.CCSize, this.CCSize, size(Points,1)));
    
%     windowKernel = parallel.gpu.CUDAKernel('windowDispROI.ptx', 'windowDispROI.cu');
%     windowKernel.ThreadBlockSize = ThreadBlockSize;
%     windowKernel.GridSize = [ceil(this.CCSize / ThreadBlockSize(1)), ceil(this.CCSize / ThreadBlockSize(2)), ceil(size(cPoints,1) / ThreadBlockSize(3))];
%     
%     [cuda_roi_1, cuda_roi_2] = feval(windowKernel, cuda_roi_1, cuda_roi_2, cPaddedTemplateGrad, cPaddedTargetGrad, this.cHann2D,...
%             cPoints, cDispX, cDispY, this.HalfCCSize, size(cPoints,1), size(cPaddedTemplateGrad,2), size(cPaddedTemplateGrad,1));
    
    for i = 1:size(Points,1)
        roi_1(:,:,i) = PaddedTemplateGrad(Points(i,2) - this.HalfCCSize:Points(i,2) + this.HalfCCSize, Points(i,1) - this.HalfCCSize:Points(i,1) + this.HalfCCSize) .* this.cHann2D;
        roi_2(:,:,i) = PaddedTargetGrad(Points(i,2) + DispY(i) - this.HalfCCSize:Points(i,2) + DispY(i) + this.HalfCCSize, Points(i,1) + DispX(i) - this.HalfCCSize:Points(i,1) + DispX(i) + this.HalfCCSize) .* this.cHann2D;
    end
        
    roi_1 = roi_1 - mean(roi_1, [1, 2]);
    roi_2 = roi_2 - mean(roi_2, [1, 2]);
    
    roi_1 = fft2(roi_1);
    roi_2 = fft2(roi_2);
    
    cuda_phase = angle(fftshift(fftshift(roi_1 .* conj(roi_2), 1), 2));
    cuda_phase = cuda_phase(3:end-3, 3:end-3, :);
    
    for i = 1:size(cuda_phase, 3)
        cuda_phase(:,:,i) = medfilt2(cuda_phase(:,:,i), [3, 3]);
    end
    
    boundary = round(this.CCSize * 0.25);
    cuda_central_phase = cuda_phase(boundary:end-boundary, boundary:end-boundary, :);
    
    [dx, dy] = fitPlane(this.cA, cuda_central_phase);
    
    dx = gather(dx);
    dy = gather(dy);
    
    dx = round(dx, 4);
    dy = round(dy, 4);
    
    subX = (this.CCSize / (2 * pi)) .* dx';
    subY = (this.CCSize / (2 * pi)) .* dy';
end

function [dx, dy] = fitPlane (A, Z)
    Zs = reshape(Z, [], size(Z,3));
    
    % Fit (ax + by + c) to the data.
    x = A \ Zs;
    
    dx = x(1,:);
    dy = x(2,:);
end