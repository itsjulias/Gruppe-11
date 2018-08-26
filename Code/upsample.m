function [img_up] = upsample(img,steps)
%UPSAMPLE Summary of this function goes here
%   Detailed explanation goes here
for i=1:steps
    [row,col]=size(img);
    if(mod(col,2) == 1 & i<=steps-2)
        col_new = col+0.5;
    else
        col_new = col;
    end
    img = interp2((1:col)',1:row,img,(1:0.5:col_new)',1:0.5:row);
    size(img) %for debugging purpose
end
img_up = img;
end

