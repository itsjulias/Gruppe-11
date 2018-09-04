function [img_segmented] = detect_homogeneous_area(img,sigma_thd,window_length,class_labels)
%DETECT_HOMOGENEOUS_AREA Summary of this function goes here
%   Detailed explanation goes here
if(mod(window_length,2)==0)
    error('window_length must be odd');
end

img = double(img);

img_segmented = zeros(size(img));

border_dist = floor(window_length/2);

for zeile = border_dist+1:(size(img,1)-border_dist-1)
    for pixel = border_dist+1:(size(img,2)-border_dist-1)
        %Bildausschnitt extrahieren
        box = img(zeile+(-border_dist:+border_dist),pixel+(-border_dist:+border_dist));
        
        std_abweichung = std(box(:));
        mean_img = mean(box(:));
        
        if std_abweichung < sigma_thd 
            % Abstand zu Klassen berechnen
            dist_class = abs(mean_img-class_labels);
            [val,idx] = min(dist_class,[],2);
            img_segmented(zeile,pixel) = class_labels(idx); 
        end
    end
end



img_segmented = uint8(img_segmented);

end

