function depth_map = depth_estimation_SH(img1_rectified,img2_rectified,K,T,offset_x_pixel,varargin)
%DEPTH_MAP Summary of this function goes here
%   Detailed explanation goes here

% % Toleranz Verschiebung der Liniensuche in y-Richtung
% tol = 0;
% Depth Map Matrix: Zeile entspricht x-Koordinate, Spalte entspricht
% y-Koordinate, Wert enspricht der Tiefe.
% Vorbelegung mit Nullen
[merkmale_img1, merkmale_img2] = ...
    extract_1D_features(img1_rectified,img2_rectified,...
    'segment_length',7,'tau',0.05e6,'min_dist',10,'N',10000,...
    'size_filter_kernel',15,'do_plot',true);
offset = 2057;
disparity_L = zeros(size(img1_rectified));
mem_kor = [];
for i = 1:size(img1_rectified,1)
    stempel1 = merkmale_img1(2,:)==i;
    stempel2 = merkmale_img2(2,:)==i;
    kor = punkt_korrespondenzen(img1_rectified,img2_rectified,...
        [merkmale_img1(1,stempel1); merkmale_img1(2,stempel1)],...
        [merkmale_img2(1,stempel2); merkmale_img2(2,stempel2)],...
        'window_length',25,'min_corr',0.5,'do_plot',false);
    if(~isempty(kor))
        disparity_L(i,kor(1,:)) = (kor(1,:)+offset_x_pixel-kor(3,:));
        mem_kor = [mem_kor, kor];
    end
end
depth_map = disparity_L;
% depth_map = zeros(size(img1_rectified,2),size(img1_rectified,1));
% for i = 25:size(img1_rectified,1)
%     leng = size(img1_rectified,2);
%     zeilenkoord = [1:leng;i*ones(1,leng)];
%     kor = punkt_korrespondenzen(img1_rectified,img2_rectified,...
%         zeilenkoord,zeilenkoord,...
%         'window_length',25,'min_corr',0.5);
%     depth = (K*T)/(kor(1:2,:)-kor(3:4,:));
%     for j = 1:size(kor,2)
%         p = kor(1:2,i);
%         depth_map(p(1),p(2)) = depth(j);
%     end
% end
% figure;
% [X,Y] = meshgrid(1:size(depth_map,1),1:size(depth_map,2));
% surf(X,Y,depth_map)
% end
%
