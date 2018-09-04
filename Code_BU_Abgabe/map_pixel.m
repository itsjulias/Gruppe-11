function [img_mapped_L,img_mapped_R] = map_pixel(img_L_rect,img_R_rect,...
    depth_map_L,depth_map_R,p,min_L_x,max_L_x,...
    min_R_x,max_R_x)
%MAP_PIXEL Summary of this function goes here
%   Detailed explanation goes here
[rows, cols] = size(img_L_rect);
img_mapped_L = zeros(rows,length(min_L_x:max_L_x));

%Pixel verschieben sich nach links
translation_L = -p;
for row=1:rows
    for pix = 1:cols
        map_pixel = round(pix+translation_L./depth_map_L(row,pix)-min_L_x);
        if(map_pixel<1 || map_pixel>max_L_x-min_L_x)
            continue;
        else
            img_mapped_L(row,map_pixel) = img_L_rect(row,pix);
        end
    end
end


[rows, cols] = size(img_R_rect);
img_mapped_R = zeros(rows,length(min_R_x:max_R_x));

% Pixel verschieben sich nach rechts
translation_R = p;
for row = 1:rows
    for pix = 1:cols
        map_pixel = round(pix+translation_R./depth_map_R(row,pix)-min_R_x);
        if(map_pixel<1 || map_pixel>max_R_x-min_R_x)
            continue;
        else
            img_mapped_R(row,map_pixel) = img_R_rect(row,pix);
        end
    end
end
img_mapped_L = uint8(img_mapped_L);
img_mapped_R = uint8(img_mapped_R);

end