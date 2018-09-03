function [disparity_map] = sad_scanline(img1,img2,window_length,offset_pix,search_r,search_l,direction)
%NCC_SCANLINE Summary of this function goes here
%   Detailed explanation goes here

% Überprüfen ob beide Bilder gleiche Anzahl an Zeilen haben
size_img1 = size(img1);
size_img2 = size(img2);

if(size_img1(1) ~=size_img2(1))
    error('Bilder haben nicht die gleiche Anzahl an Zeilen')
end

if(mod(window_length,2)==0)
    error('window_length must be odd')
end

if(~strcmp(direction,'RL')&&~strcmp(direction,'LR'))
    error('direction must be either RL or LR')
end

if(search_r<0 || search_l<0)
    error('search_r oder search_l must be greater than 0')
end

% if(size_img1(2)<size_img2(2))
%     img2 = img2(1:size_img1(1),1:size_img1(2));
% elseif(size_img1(2)>size_img2(2))
%     img1 = img1(1:size_img2(1),1:size_img2(2));
% end

h = waitbar(0,['Disparity map generation from ', direction(1) ' to ' direction(2)]);
set(h,'Name','Disparity generation progress');

img1 = double(img1);
img2 = double(img2);
size_x = size(img1,2);
size_y = size(img1,1);

size_x2 = size(img2,2);
size_y2 = size(img2,1);

% Rand für Bildsegmentauswahl zur Verfügung stellen
border_dist = floor(window_length/2);

% Mpt1-Vektor (Merkmalspunkte im Bild 1 entlang x-Achse)
Mpt1 = (border_dist+1):(size_x-border_dist-1);

% Mpt2-Vektor (Merkmalspunkte im Bild 2 entlang x-Achse)
Mpt2 = (border_dist+1):(size_x2-border_dist-1);
Mpt2_off = Mpt2-offset_pix;

% Für Schleife benötigte Parameter
window_length_squared = window_length^2;

% Bestimme Disparity-Band
if(strcmp(direction,'LR'))
    max_disparity =  search_r; %abs(max(max(rough_disparity(low_y:up_y,:)))+offset_pix);
    min_disparity = search_l; %abs(min(min(rough_disparity(low_y:up_y,:)))+offset_pix);
    disparity_map = zeros(size_y,size_x);

elseif(strcmp(direction,'RL'))
    max_disparity =  search_l; %abs(max(max(rough_disparity(low_y:up_y,:)))+offset_pix);
    min_disparity = search_r; %abs(min(min(rough_disparity(low_y:up_y,:)))+offset_pix);
    disparity_map = zeros(size_y2,size_x2);

end

% Zeilenweise wird SAD berechnet
for zeile = (border_dist+1):1:(size_y-border_dist-1)
    % Berechnen der oberen und unteren Fenster indizes pro Zeile
    low_y = zeile-border_dist;
    up_y = zeile+border_dist;
    low_x = Mpt1-border_dist;
    up_x = Mpt1+border_dist;
    
    low2_y = zeile-border_dist;
    up2_y = zeile+border_dist;
    low2_x = Mpt2-border_dist;
    up2_x = Mpt2+border_dist;
    
    
    % Über jedes Pixel der Zeile von Bild1 wird Box extrahiert
    for i = 1:length(Mpt1)
        % Bildausschnitt extrahieren
        boxes_img1(:,i) = (reshape(img1(low_y:up_y,low_x(i):up_x(i)),[window_length_squared,1]));...
%             -mean(mean(img1(low_y:up_y,low_x(i):up_x(i))));
    end
    
    % Über jedes Pixel der Zeile von Bild1 wird Box extrahiert
    for i = 1:length(Mpt2)
        % Bildausschnitt extrahieren
        boxes_img2(:,i) = (reshape(img2(low2_y:up2_y,low2_x(i):up2_x(i)),[window_length_squared,1]));...
%             -mean(mean(img2(low_y:up_y,low_x(i):up_x(i))));
    end
    
    if(strcmp(direction,'LR'))
        % Für jedes Pixel im linken Bild wird entsprechende Korrespondenz
        % im rechten Bild gesucht und Disparitaet dazu ermittelt
        for i = 1:1:length(Mpt1)
            SAD = sum(abs(boxes_img1(:,i)...
                -boxes_img2(:,max(1,i-min_disparity):min(length(Mpt2),i+max_disparity)))); %Suche in eingeschränktem Bereich (max() und min() damit nicht über validen Bereich hinausgesucht wird)
            %
            % SAD = sum(abs(boxes_img1(:,i)-boxes_img2)); %Suche durch komplettes Bild  
            % SAD = sum((boxes_img1(:,i)-boxes_img2).^2); %Suche mit SSD durch komplettes Bild
            [val, idx] = sort(SAD,'ascend');
            idx_min = min(idx(1)+max(1,i-min_disparity),length(Mpt2_off));%ceil(mean(idx(val == min(val))));
            if(~isempty(idx))
                new_disparity =max(Mpt1(i)- Mpt2_off(idx_min),0); 
                disparity_map(zeile,Mpt1(i))= new_disparity;
            end
        end
    elseif(strcmp(direction,'RL'))
        % Für jedes Pixel im rechten Bild wird entsprechende Korrespondenz
        % im linken Bild gesucht und Disparitaet dazu ermittelt
        for i = 1:1:length(Mpt2)
            SAD = sum(abs(boxes_img2(:,i)...
                -boxes_img1(:,max(1,i-min_disparity):min(length(Mpt1),i+max_disparity)))); %Suche in eingeschränktem Bereich (max() und min() damit nicht über validen Bereich hinausgesucht wird)
            % SAD = sum(abs(boxes_img1(:,i)-boxes_img2)); %Suche durch komplettes Bild  
            % SAD = sum((boxes_img1(:,i)-boxes_img2).^2); %Suche mit SSD durch komplettes Bild
            [val, idx] = sort(SAD,'ascend');
            idx_min = min(idx(1)+max(1,i-min_disparity),length(Mpt1));%ceil(mean(idx(val == min(val))));
            if(~isempty(idx))
                new_disparity =max(Mpt1(idx_min)- Mpt2_off(i),0);
                disparity_map(zeile,Mpt2(i))= new_disparity;
            end
        end
    end
    waitbar(zeile/(size_y-border_dist-1));
end
close(h)
end

