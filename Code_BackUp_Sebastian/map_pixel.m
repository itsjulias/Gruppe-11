function [img_mapped_L,img_mapped_R,img_mapped] = map_pixel(img_L_rect,img_R_rect,...
    depth_map_L,depth_map_R,p,min_L_x,max_L_x,...
    min_R_x,max_R_x)
%MAP_PIXEL Summary of this function goes here
%   Detailed explanation goes here
[rows, cols] = size(img_L_rect);
img_mapped_L = zeros(rows,length(min_L_x:max_L_x));

%% Mapping der Pixel aus dem linken Bild
% Beim mapping der Pixel aus dem linken Bild in die virtuellen Ansicht
% verschieben sich die Pixel im Bild nach links. Deshalb muss für den
% Translationsvektor -p in den rektifizierten Kooridnaten angenommen
% werden.

translation_L = -p; % Pixel verschieben sich nach links
for row=1:rows % Iteriere durch jede Zeile
    for pix = 1:cols 
        % Berechne für jedes Pixel das Zielpixel in der virtuellen Ansicht
        map_pixel = round(pix+translation_L./depth_map_L(row,pix)-min_L_x);
        if(map_pixel<1 || map_pixel>max_L_x-min_L_x)
            % Wenn das map_pixel (Zielpixel) außerhalb der Grenzen fuer die
            % virtuelle Ansicht liegt -> Fortfahren ohne Eintrag.
            continue;
        else 
            % Wenn das map_pixel innerhalb der Grenzen fuer die virtuelle
            % fuehre mapping durch.
            img_mapped_L(row,map_pixel) = img_L_rect(row,pix);
        end
    end
end


[rows, cols] = size(img_R_rect);
img_mapped_R = zeros(rows,length(min_R_x:max_R_x));
%% Mapping der Pixel aus dem rechten Bild
% Beim mapping der Pixel aus dem rechten Bild in die virtuellen Ansicht
% verschieben sich die Pixel im Bild nach rechts. Deshalb muss für den
% Translationsvektor 1-p (p>0) in den rektifizierten Kooridnaten angenommen
% werden.

translation_R = 1-p;% Pixel verschieben sich nach rechts
for row = 1:rows
    for pix = 1:cols
        % Berechne für jedes Pixel das Zielpixel in der virtuellen Ansicht
        map_pixel = round(pix+translation_R./depth_map_R(row,pix)-min_R_x);
        if(map_pixel<1 || map_pixel>max_R_x-min_R_x)
            % Wenn das map_pixel (Zielpixel) außerhalb der Grenzen fuer die
            % virtuelle Ansicht liegt -> Fortfahren ohne Eintrag.
            continue;
        else
            % Wenn das map_pixel innerhalb der Grenzen fuer die virtuelle
            % fuehre mapping durch.
            img_mapped_R(row,map_pixel) = img_R_rect(row,pix);
        end
    end
end



%% Mapping der Pixel aus linkem & rechten Bild in virtuelle Ansicht
% Die Translationen wurden zuvor (s. vorherige Abschnitte) festgelegt.
translation_L = -p; % Pixel verschieben sich nach links
translation_R = 1-p;% Pixel verschieben sich nach rechts

% Da Rueckprojektion immer vom linken Bild ausgehend stattfindet, wird die
% Groesse des linken Bilds verwendet um spätere
% size-Ungenauigkeiten/Unterschiede bei der Interpolation zu vermeiden
rows = size(img_L_rect,1);%min(size(img_R_rect,1),size(img_L_rect,1));
% Für die for-Schleife wird jedoch das kleinste Bild beruecksichtigt, das
% sonst der Fehler "Index exceeds matrix dimensions." auftritt
cols = min(size(img_R_rect,2),size(img_L_rect,2));%size(img_L_rect,2);

img_mapped = zeros(rows,length(min_L_x:max_L_x));
for row = 1:rows
    for pix = 1:cols
        % Selbes Mapping wie in obigen Abschnitten, nur dass linkes Bild
        % vor rechtem Bild gemapped wird.
        map_pixel_L = round(pix+translation_L./depth_map_L(row,pix)-min_L_x);
        if(map_pixel_L<1 || map_pixel_L>size(img_mapped,2))%max_L_x-min_L_x)
            continue;
        else
            img_mapped_L(row,map_pixel_L) = img_L_rect(row,pix);
        end
        map_pixel_R = round(pix+translation_R./depth_map_R(row,pix)-min_R_x);
        if(map_pixel_R<1 || map_pixel_R>size(img_mapped,2) || img_mapped(row,map_pixel_R)~=0)
            continue;
        else
            img_mapped(row,map_pixel_R) = img_R_rect(row,pix);
        end
    end
end

% Gebe die einzeln gemappten Bilder zurueck
img_mapped_L = uint8(img_mapped_L);
img_mapped_R = uint8(img_mapped_R);
img_mapped = uint8(img_mapped);

end