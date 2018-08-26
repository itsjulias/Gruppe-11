function [depth_map,depth_map2,disparity_map] = depth_estimation(img1_rectified,img2_rectified,K,T,offset_x_pixel,varargin)
%DEPTH_MAP Summary of this function goes here
%   Detailed explanation goes here
disparity_map = zeros(min(size(img1_rectified),size(img2_rectified)));
for i = 3:-1:1
    %downsampling of images
    [img1_dwn, img2_dwn] = downsample(img1_rectified,img2_rectified,i);
    %NCC disparity map generation
    disparity_map_dwn = ncc_scanline(img1_dwn,img2_dwn,15,offset_x_pixel/(2^i),0.8);
    disparity_map_up(:,:,i) = upsample((2^i)*disparity_map_dwn,i);
    figure
    surf(disparity_map_up(:,:,i),'LineStyle','none')
    view(0,-90);
    disparity_map = max(disparity_map,disparity_map_up(:,:,i));
end
% disparity_map = sum(disparity_map_up(:,:,:),3)/size(disparity_map_up,3);
depth_map2 = 1./max(disparity_map,zeros(size(disparity_map)));
depth_map = norm(K*T)./max(disparity_map,zeros(size(disparity_map)));
figure
surf(depth_map,'LineStyle','none')
view(0,-90);

end