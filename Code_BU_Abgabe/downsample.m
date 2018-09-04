function [img1_dwn,odd_even_xy] = downsample(img1,steps,get_every_layer)
%DOWNSAMPLE Summary of this function goes here
%   Detailed explanation goes here
get_every_layer = nargin == 4 && strcmp(get_every_layer,'get_every_layer');

size_filter_kernel = 5;
sigma = floor(size_filter_kernel/2)/3;
w = 1/(sqrt(2*pi)*sigma)*exp(-(-floor(size_filter_kernel/2):1:floor(size_filter_kernel/2)).^2/(2*sigma^2));
img1_filt = double(img1);
odd_even_xy = [];
for i=1:steps
    odd_even_xy = [odd_even_xy; mod(size(img1_filt),2)]; % for upsampling 
    size(img1_filt)
    img1_filt(ceil(size_filter_kernel/2):(end-floor(size_filter_kernel/2)),...
              ceil(size_filter_kernel/2):(end-floor(size_filter_kernel/2))) = ...
        conv2(img1_filt,w'*w,'valid'); 
    img1_filt = img1_filt(1:2:end,1:2:end);
    if get_every_layer
        img1_dwn{i} = uint8(img1_filt);
    end
end
img1_dwn = uint8(img1_filt);

if ~get_every_layer
    img1_dwn = uint8(img1_filt);
end

end
