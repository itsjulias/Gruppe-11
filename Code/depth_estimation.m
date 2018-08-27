function depth_map = ...
    depth_estimation(img1_rectified,img2_rectified,K,T,offset_x_pixel,...
    d_cut_up,d_cut_down,varargin)
%DEPTH_MAP Summary of this function goes here
%   Detailed explanation goes here
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
end

dist_side = 5;
% Depth map nach oben erweitern mit der bisherigen obersten Zeile
depth_map = [repmat(depth_map(dist_side+1,:),d_cut_up,1); depth_map];
% Depth map nach unten erweitern mit der bisherigen untersten Zeile
depth_map = [depth_map;repmat(depth_map(end-dist_side,:),d_cut_down,1)];
end
