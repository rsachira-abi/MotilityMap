function [intX, intY, intError] = calcIntegerShift (this, cPaddedTemplateGrad, cPaddedTargetGrad, cPoints, ThreadBlockSize)
    %keyboard;
    cuda_roi_1 = complex(zeros(this.CCSize, this.CCSize, size(cPoints,1), 'gpuArray'));
    cuda_roi_2 = complex(zeros(this.CCSize, this.CCSize, size(cPoints,1), 'gpuArray'));
    
    windowKernel = parallel.gpu.CUDAKernel('windowROI.ptx', 'windowROI.cu');
    windowKernel.ThreadBlockSize = ThreadBlockSize;
    windowKernel.GridSize = [ceil(this.CCSize / ThreadBlockSize(1)), ceil(this.CCSize / ThreadBlockSize(2)), ceil(size(cPoints,1) / ThreadBlockSize(3))];
    
    [cuda_roi_1, cuda_roi_2] = feval(windowKernel, cuda_roi_1, cuda_roi_2, cPaddedTemplateGrad, cPaddedTargetGrad, this.cHamming2D,...
            cPoints, this.HalfCCSize, size(cPoints,1), size(cPaddedTemplateGrad,2), size(cPaddedTemplateGrad,1));
    
    cuda_roi_1 = cuda_roi_1 - mean(cuda_roi_1, [1, 2]);
    cuda_roi_2 = cuda_roi_2 - mean(cuda_roi_2, [1, 2]);
    
    cuda_roi_1 = fft2(cuda_roi_1);
    cuda_roi_2 = conj(fft2(cuda_roi_2));
    
    CC = abs(fftshift(fftshift(ifft2(cuda_roi_1 .* cuda_roi_2), 1), 2));
    
    maxCC = max(CC, [], [1, 2]);

    % Avoid NaNs. Think of a better solution later.
    % Fill all NaN rois with 1s so that they are omited due to high
    % error in the later steps
    zero_max_logical = (maxCC == 0);
    CC(:,:,zero_max_logical) = 1;
    maxCC(zero_max_logical) = 1;

    CC = CC ./ maxCC;
    
    ind = find(CC == 1);
    [y, x, z] = ind2sub(size(CC), ind);
    
    x = gather(x);
    y = gather(y);
    z = gather(z);

    % Take the mean of the x y values of 1s as the (x, y) of that z
    y = round(splitapply(@mean, y, z));
    x = round(splitapply(@mean, x, z));

    x = gpuArray(x);
    y = gpuArray(y);
%     [~, ia, ~] = unique(z);
%     x = x(ia);
%     y = y(ia);
    
    intX = ((size(CC,2) - 1) / 2) - x + 1;
    intY = ((size(CC,1) - 1) / 2) - y + 1;

    AroundPeakKernel = parallel.gpu.CUDAKernel('around_peak.ptx', 'around_peak.cu');
    AroundPeakKernel.ThreadBlockSize = windowKernel.ThreadBlockSize;
    AroundPeakKernel.GridSize = windowKernel.GridSize;
    
    CC = feval(AroundPeakKernel, CC, int32(x), int32(y), this.SquareSize, this.PeakThreshold, this.HalfCCSize, size(cPoints,1));
    intError = squeeze(sum(CC, [1, 2]));
end