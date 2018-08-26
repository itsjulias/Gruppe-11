function [depth_map,depth_map2,disparity_map] = depth_estimation(img1_rectified,img2_rectified,K,T,offset_x_pixel,varargin)
%DEPTH_MAP Summary of this function goes here
%   Detailed explanation goes here
version =2;
if version ==1
    disparity_map = zeros(min(size(img1_rectified),size(img2_rectified)));
    for i = 2:-1:1
        %downsampling of images
        [img1_dwn, img2_dwn,odd_even_xy] = downsample(img1_rectified,img2_rectified,i);
        %NCC disparity map generation
        disparity_map_dwn = sad_scanline(img1_dwn,img2_dwn,15,offset_x_pixel/(2^i),0.75);
        disparity_map_up(:,:,i) = upsample((2^i)*disparity_map_dwn,i,odd_even_xy);
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
elseif(version==2)
    for i = 2:-1:2
        %downsampling of images
        [img1_dwn, img2_dwn,odd_even_xy] = downsample(img1_rectified,img2_rectified,i);
        %NCC disparity map generation
        depth_map_dwn = depth_estimation4(img1_dwn,img2_dwn, offset_x_pixel/(2^i));
        depth_map_up(:,:,i) = upsample(1/(2^i)*depth_map_dwn,i,odd_even_xy);
        figure
        surf(depth_map_up(:,:,i),'LineStyle','none')
        view(0,-90);
        depth_map = depth_map_up(:,:,i);%min(depth_map_up,depth_map_up(:,:,i));
        depth_map2 = 0;
        disparity_map = 0;
        
    end
end