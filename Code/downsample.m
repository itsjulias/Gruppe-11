function [img1_dwn,img2_dwn] = downsample(img1,img2,steps)
%DOWNSAMPLE Summary of this function goes here
%   Detailed explanation goes here
size_filter_kernel = 5;
sigma = floor(size_filter_kernel/2)/4;
w = 1/(sqrt(2*pi)*sigma)*exp(-(-floor(size_filter_kernel/2):1:floor(size_filter_kernel/2)).^2/(2*sigma^2));
img1_filt = double(img1);
img2_filt = double(img2);

for i=1:steps
    img1_filt = conv2(img1_filt,w'*w,'same');
    img2_filt = conv2(img2_filt,w'*w,'same');
    img1_filt = img1_filt(1:2:end,1:2:end);
    img2_filt = img2_filt(1:2:end,1:2:end);   
    size(img2_filt) % for debugging purpose
end
img1_dwn = uint8(img1_filt);
img2_dwn = uint8(img2_filt);

end

