close all;
clear all;
load('zwspeicher');

disparity_L = depth_estimation_SH(img1_rectified,img2_rectified,K,T);
figure
imshow(disparity_L,jet)

disparity_L_vL = disparity_L;
for i= 1:size(disparity_L_vL,1)
    for j=2:size(disparity_L_vL,2)-1
        if(disparity_L_vL(i,j)==0)
            disparity_L_vL(i,j) = disparity_L_vL(i,j-1);
        else
            continue;
        end
    end
end

disparity_L_vR = disparity_L;
for i= 1:size(disparity_L_vR,1)
    j=size(disparity_L_vR,2)-1;
    while j~=2
        if(disparity_L_vR(i,j)==0)
            disparity_L_vR(i,j) = disparity_L_vR(i,j+1);
            j = j-1;
        else
            j = j-1;
        end
    end
end

% disparity_L_interp = disparity_L;
% size_disp = size(disparity_L_interp);
% for i= 1:size_disp(2)
%     [val, idx] = find(disparity_L_interp(:,i) ~= 0);
%     if(size(idx,1)>2)
%         disparity_L_interp(:,i) = interp1(idx,disparity_L_interp(idx,i),1:size_disp(1)','linear',NaN);
%     end
%     
% end

size_filter_kernel = 50;
sigma = floor(size_filter_kernel/2)/3;
w = 1/(sqrt(2*pi)*sigma)*exp(-(-floor(size_filter_kernel/2):1:floor(size_filter_kernel/2)).^2/(2*sigma^2));
surf(conv2(min(disparity_L_vR,disparity_L_vL),w'*w,'same'),'LineStyle','none')
view(0,-90);



figure
imshow(disparity_L,jet)