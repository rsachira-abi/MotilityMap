
function grad_image = grad (this, image)
    Gx = conv2(image, this.GradKernel, 'same');
    Gy = conv2(image, this.GradKernel', 'same');
    
    grad_image = Gx + Gy*1i;
end