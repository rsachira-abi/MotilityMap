
function [intX, intY, intError] = calcIntegerShift (this, PaddedTemplateGrad, PaddedTargetGrad, Points)
    %keyboard;
    roi_1 = complex(zeros(this.CCSize, this.CCSize, size(Points,1), 'gpuArray'));
    roi_2 = complex(zeros(this.CCSize, this.CCSize, size(Points,1), 'gpuArray'));
    
%     windowKernel = parallel.gpu.CUDAKernel('windowROI.ptx', 'windowROI.cu');
%     windowKernel.ThreadBlockSize = ThreadBlockSize;
%     windowKernel.GridSize = [ceil(this.CCSize / ThreadBlockSize(1)), ceil(this.CCSize / ThreadBlockSize(2)), ceil(size(cPoints,1) / ThreadBlockSize(3))];
%     
%     [cuda_roi_1, cuda_roi_2] = feval(windowKernel, cuda_roi_1, cuda_roi_2, cPaddedTemplateGrad, cPaddedTargetGrad, this.cHamming2D,...
%             cPoints, this.HalfCCSize, size(cPoints,1), size(cPaddedTemplateGrad,2), size(cPaddedTemplateGrad,1));

    for i = 1:size(Points,1)
        roi_1(:,:,i) = PaddedTemplateGrad(Points(i,2) - this.HalfCCSize:Points(i,2) + this.HalfCCSize, Points(i,1) - this.HalfCCSize:Points(i,1) + this.HalfCCSize) .* this.cHamming2D;
        roi_2(:,:,i) = PaddedTargetGrad(Points(i,2) - this.HalfCCSize:Points(i,2) + this.HalfCCSize, Points(i,1) - this.HalfCCSize:Points(i,1) + this.HalfCCSize) .* this.cHamming2D;
    end
    
    roi_1 = roi_1 - mean(roi_1, [1, 2]);
    roi_2 = roi_2 - mean(roi_2, [1, 2]);
    
    roi_1 = fft2(roi_1);
    roi_2 = conj(fft2(roi_2));
    
    CC = abs(fftshift(fftshift(ifft2(roi_1 .* roi_2), 1), 2));
    
    maxCC = max(CC, [], [1, 2]);
    CC = CC ./ maxCC;
    
    ind = find(CC == 1);
    [y, x, z] = ind2sub(size(CC), ind);
    [~, ia, ~] = unique(z);
    x = x(ia);
    y = y(ia);
    
    intX = int32(((size(CC,2) - 1) / 2) - x + 1);
    intY = int32(((size(CC,1) - 1) / 2) - y + 1);

%     AroundPeakKernel = parallel.gpu.CUDAKernel('around_peak.ptx', 'around_peak.cu');
%     AroundPeakKernel.ThreadBlockSize = windowKernel.ThreadBlockSize;
%     AroundPeakKernel.GridSize = windowKernel.GridSize;
%     
%     CC = feval(AroundPeakKernel, CC, int32(x), int32(y), this.SquareSize, this.PeakThreshold, this.HalfCCSize, size(Points,1));
    
    i_min = max([intY - this.SquareSize, ones(size(intY))], [], 2);
    i_max = min([intY + this.SquareSize, ones(size(intY)) .* size(CC,1)], [], 2);
    j_min = max([intX - this.SquareSize, ones(size(intX))], [], 2);
    j_max = min([intX + this.SquareSize, ones(size(intX)) .* size(CC,2)], [], 2);
    
    for k = 1:size(CC,3)
        a = CC(i_min(k):i_max(k), j_min(k):j_max(k), k) > this.PeakThreshold;
        intError(k) = sum(a(:));
    end
end