function [disp_map_filt] = remove_outliers(disp_map,window_length,type)
%REMOVE_OUTLIERS Summary of this function goes here
%   Detailed explanation goes here
if(mod(window_length,2) == 0)
    error('window_length must be odd')
end

[size_row,size_col] = size(disp_map);
disp_map_filt = disp_map;
border_dist = floor(window_length/2);


for row = (border_dist+1):(size_row-(border_dist+1))
    for col = (border_dist+1):(size_col-(border_dist+1))
        if(disp_map(row,col) ~= 0 && ~isnan(disp_map(row,col)) && disp_map(row,col)~= Inf)
            % Extrahiere Fenster um Pixel
            window = disp_map(row+(-border_dist:border_dist),...
                col+(-border_dist:border_dist));
            if(~isempty(window(window~=0 & ~isnan(window)& isfinite(window))))
                if(strcmp(type,'median'))
                    disp_map_filt(row,col) = median(window(window~=0 & ~isnan(window)& window~=Inf));
                elseif(strcmp(type,'mean'))
                    disp_map_filt(row,col) = mean(window(window~=0));
                end
            end
        else
            disp_map_filt(row,col) = 0;
        end
    end
end

end

