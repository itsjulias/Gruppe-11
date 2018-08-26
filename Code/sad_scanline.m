function [disparity_map] = sad_scanline(img1,img2,window_length,offset_pix,rough_disparity)
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

if(size_img1(2)<size_img2(2))
    img2 = img2(1:size_img1(1),1:size_img1(2));
elseif(size_img1(2)>size_img2(2))
    img1 = img1(1:size_img2(1),1:size_img2(2));
end

size_x = size(img1,2);
size_y = size(img1,1);

% Rand für Bildsegmentauswahl zur Verfügung stellen
border_dist = floor(window_length/2);

% Mpt1-Vektor (Merkmalspunkte im Bild 1 entlang x-Achse)
Mpt1 = (border_dist+1):(size_x-border_dist-1);

% Mpt2-Vektor (Merkmalspunkte im Bild 2 entlang x-Achse)
Mpt2 = (border_dist+1):(size_x-border_dist-1);

% Für Schleife benötigte Parameter
window_length_squared = window_length^2;
disparity_map = zeros(size_y,size_x);
disparity_map_L = zeros(size_y,size_x);
disparity_map_R = zeros(size_y,size_x);

% Zeilenweise wird SAD berechnet
for zeile = (border_dist+1):1:(size_y-border_dist-1)
    % Berechnen der oberen und unteren Fenster indizes
    low_y = zeile-border_dist;
    up_y = zeile+border_dist;
    low_x = Mpt1-border_dist;
    up_x = Mpt1+border_dist;
    
    low2_y = zeile-border_dist;
    up2_y = zeile+border_dist;
    low2_x = Mpt1-border_dist;
    up2_x = Mpt1+border_dist;
    
    % Bestimme Disparity-Band
    max_disparity =  1000; %abs(max(max(rough_disparity(low_y:up_y,:)))+offset_pix);
    min_disparity = -300; %abs(min(min(rough_disparity(low_y:up_y,:)))+offset_pix);
    
    % Über jedes Pixel der Zeile von Bild1 wird Box extrahiert
    for i = 1:length(Mpt1)
        % Bildausschnitt extrahieren
        boxes_img1(:,i) = reshape(img1(low_y:up_y,low_x(i):up_x(i)),[window_length_squared,1])...
            -mean(mean(img1(low_y:up_y,low_x(i):up_x(i))));
    end
    
    % Über jedes Pixel der Zeile von Bild1 wird Box extrahiert
    for i = 1:length(Mpt2)
        % Bildausschnitt extrahieren
        boxes_img2(:,i) = reshape(img2(low_y:up_y,low_x(i):up_x(i)),[window_length_squared,1])...
            -mean(mean(img2(low_y:up_y,low_x(i):up_x(i))));
    end
    
%     pix_R_old = 1;
    pix_R = 1;
    cnt = 1;
    for i = 1:1:length(Mpt1)
%         SAD = sum(abs(boxes_img1(:,i)...
                 -boxes_img2(:,max(1,i-min_disparity):min(length(Mpt2),i+max_disparity))));
%       
        SAD = sum(abs(boxes_img1(:,i)-boxes_img2));
%         SAD = sum((boxes_img1(:,i)-boxes_img2).^2);
        [val idx] = sort(SAD,'ascend');
        idx_min = ceil(mean(idx(val == min(val))));
        if(~isempty(idx))
%             pix_R_old = pix_R;
%             pix_R = idx_min(1);%idx_min(ceil(rand*length(idx_min)));
            
            new_disparity =max(min(Mpt1(i)- Mpt2(idx_min)+offset_pix,500),50);
           
%             pix_R+offset_pix-(i+border_dist);
%             if(abs(pix_R_old-pix_R)<10)
%             if((abs(new_disparity-disparity_map_L(zeile,pix-1))<min_disparity) & (pix ~=1))
%                 disparity_map_L(zeile,i) = disparity_map_L(zeile,i-1);
%             else
            disparity_map_L(zeile,i)= new_disparity;%-size_x;
%             end
        end
        cnt = cnt+1;
    end
    
    %     for pix = 1:2:length(Mpt2)
    % %         SAD = sum(abs(boxes_img2(:,pix)...
    % %                       -boxes_img1(:,max(1,pix-min_disparity):min(length(Mpt2),pix+max_disparity))));
    %         SAD = sum(abs(boxes_img2(:,pix)-boxes_img1));
    %
    %         [val idx] = sort(SAD,'ascend');
    %         idx_min = idx(val == min(val));
    %          if(~isempty(idx))
    %              pix_disparity = idx_min(1);%idx_min(ceil(rand*length(idx_min)));
    %             disparity_map_R(zeile,pix) = -(pix+border_dist-pix_disparity);
    %         end
    %     end
    
   disp(zeile);
end
disparity_map = disparity_map_L;%0.5*(disparity_map_R+disparity_map_L);
end

% %% Merkmalsvorbereitung
% [length_y_im1, length_x_im1] = size(Im1);
% [length_y_im2, length_x_im2] = size(Im2);
%
% % Alle Merkmalspunkte, die zu nahe am Rand liegen werden nicht
% % berücksichtigt/entfernt
% Mpt1 = Mpt1(:,(Mpt1(1,:)+floor(window_length/2)<=length_x_im1) & (Mpt1(1,:)-floor(window_length/2)>=1)...
%     & (Mpt1(2,:)+floor(window_length/2)<=length_y_im1) & (Mpt1(2,:)-floor(window_length/2)>=1));
% Mpt2 = Mpt2(:,(Mpt2(1,:)+floor(window_length/2)<=length_x_im2) & (Mpt2(1,:)-floor(window_length/2)>=1)...
%     & (Mpt2(2,:)+floor(window_length/2)<=length_y_im2) & (Mpt2(2,:)-floor(window_length/2)>=1));
%
% %% Merkmalsvorbereitung
% [length_y_im1, length_x_im1] = size(Im1);
% [length_y_im2, length_x_im2] = size(Im2);
%
% % Alle Merkmalspunkte, die zu nahe am Rand liegen werden nicht
% % berücksichtigt/entfernt
% Mpt1 = Mpt1(:,(Mpt1(1,:)+floor(window_length/2)<=length_x_im1) & (Mpt1(1,:)-floor(window_length/2)>=1)...
%     & (Mpt1(2,:)+floor(window_length/2)<=length_y_im1) & (Mpt1(2,:)-floor(window_length/2)>=1));
% Mpt2 = Mpt2(:,(Mpt2(1,:)+floor(window_length/2)<=length_x_im2) & (Mpt2(1,:)-floor(window_length/2)>=1)...
%     & (Mpt2(2,:)+floor(window_length/2)<=length_y_im2) & (Mpt2(2,:)-floor(window_length/2)>=1));
%
% if(isempty(Mpt1))
%     Korrespondenzen = [];
%     return
% end
%
% % Einschränken auf voraussichtlichen Disparitätsbereich
% Mpt2 = Mpt2(:,(Mpt2(1,:)<max(Mpt1(1,:)+max_disparity)) &...
%     (Mpt2(1,:)>min(Mpt1(1,:)+min_disparity)));
%
% %% Normierung
% section = zeros(window_length^2,length(Mpt1(1,:)));
% section2 = zeros(window_length^2,length(Mpt2(1,:)));
% windows2 = zeros(window_length,window_length,length(Mpt2(1,:)));
%
% low_y = (Mpt1(2,:)-window_dist);
% up_y = (Mpt1(2,:)+window_dist);
% low_x = (Mpt1(1,:)-window_dist);
% up_x = (Mpt1(1,:)+window_dist);
%
% low2_y = (Mpt2(2,:)-window_dist);
% up2_y = (Mpt2(2,:)+window_dist);
% low2_x = (Mpt2(1,:)-window_dist);
% up2_x = (Mpt2(1,:)+window_dist);
%
%
%
% for i=1:length(Mpt2(1,:))
%     % Bildausschnitt extrahieren
%     windows2(:,:,i) = Im2(low2_y(i):up2_y(i),low2_x(i):up2_x(i));
% end
%
% section2 = reshape(windows2,[window_length_squared,length(Mpt2(1,:))]);
%
%
% %% SAD Brechnung
% SAD = sum(abs(section-section2));
% [val, sorted_index] = sort(SAD,'ascend');
%
%
% %% Korrespondenz
% Korrespondenzen = [Mpt1; sorted_index(1)+window_dist; Mpt1(end,:)];
%
% if(do_plot)
%     figure;
%     im1 = imshow(I1); hold all;
%     im2 = imshow(I2);
%     im1.AlphaData = 0.5;
%     im2.AlphaData = 0.5;
%
%     for i=1:length(Korrespondenzen(1,:))
%         plot(Korrespondenzen(1,i), Korrespondenzen(2,i),'*r')
%         plot(Korrespondenzen(3,i), Korrespondenzen(4,i),'*g')
%         plot([Korrespondenzen(1,i), Korrespondenzen(3,i)], ...
%             [Korrespondenzen(2,i), Korrespondenzen(4,i)],'-b')
%     end
%
% end
%
% end
%
