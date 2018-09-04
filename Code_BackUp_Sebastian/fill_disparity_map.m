function [disparity_map_new] = fill_disparity_map(disparity_map,stamp,window_length,corner_dist)
%FILL_DISPARITY_MAP fills disparity map within stamp
% Entferne Randeffekte in Disparity Map und Stempel


disparity_map_new =zeros(size(disparity_map));
stamp_new = zeros(size(stamp));
disparity_map_new(corner_dist:(end-corner_dist),corner_dist:(end-corner_dist)) = ...
    disparity_map(corner_dist:(end-corner_dist),corner_dist:(end-corner_dist));
stamp_new(corner_dist:(end-corner_dist),corner_dist:(end-corner_dist)) =...
    stamp(corner_dist:(end-corner_dist),corner_dist:(end-corner_dist));

border_dist = floor(window_length/2);
finish_flag = (sum(sum(stamp_new(border_dist+1:end-border_dist,...
    border_dist+1:end-border_dist))) == 0);
cnt = 0;
while ~finish_flag
    %Solange NULLEN im Bild/Disparity-Map vorhanden sind fuelle diese
    %Zeilenweise und von links nach rechts auf mittels Median-Filter
    for zeile = 1:size(disparity_map,1)%(border_dist+1):(size(disparity_map,1)-border_dist)
        for pix = 1:size(disparity_map,2)% (border_dist+1):(size(disparity_map,2)-border_dist)
            
            if stamp_new(zeile,pix) ~=0
                % Stempel gibt Bereich vor, fuer welchen Nullen durch Werte
                % ersetzt werden dürfen
                
                %Window ist quadratischer Bildausschnitt, am Bildand wird
                %der jeweils ins Nichts überlappende Bereich nicht
                %ins Fenster uebertragen (z.B. bei 5x5-Fenster beim Pixel
                %(1,1) erhält man nur ein 3x3- Fenster, beim Pixel (100,1)
                %nur ein 5x3-Fenster usw.
                zeilen_ausschnitt = (-border_dist:border_dist)+zeile;
                zeilen_ausschnitt = zeilen_ausschnitt(zeilen_ausschnitt>0 & zeilen_ausschnitt<size(disparity_map,1));
                pixel_ausschnitt = (-border_dist:border_dist)+pix;
                pixel_ausschnitt = pixel_ausschnitt(pixel_ausschnitt>0 & pixel_ausschnitt<size(disparity_map,2));
                window =  disparity_map_new(zeilen_ausschnitt,pixel_ausschnitt);
                if(sum(window(:)~=0)>3)
                    % Nur wenn Fenster mehr als 3 valide Werte enthält,
                    % wird Median gebildet
                    disparity_map_new(zeile,pix) = median(window(window~=0));
                    stamp_new(zeile,pix) = 0;
                end
            end
        end
    end
    cnt = cnt+1;
    finish_flag = (sum(sum(stamp_new(border_dist+1:end-border_dist,...
        border_dist+1:end-border_dist))) == 0) || cnt == 100;
end
end

