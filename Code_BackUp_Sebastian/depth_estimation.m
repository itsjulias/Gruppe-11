function [depth_map_L, depth_map_R] = ...
    depth_estimation(img1_rectified,img2_rectified,K,T,offset_x_pixel,...
    d_cut_up_L,d_cut_down_L,d_cut_up_R,d_cut_down_R,...
    min_disparity,max_disparity,varargin)
%DEPTH_MAP 
    window_length = 5;
    % Minimale und maximale Disparität wurde über robuste Korrspondenzen
    % ermittelt. dist_safety legt fest, um wie viele Pixel die tatsächliche
    % min./max. Disparität von den ermittelten Werten abweichen darf
    dist_safety = 0.5*min(abs(min_disparity),abs(max_disparity));
    % Suchintervall festlegeben. Beispiel: Bei interv_search_left = -100
    % und interv_search_right = 200 wird 100 Pixel nach links und 200 nach
    % rechts gesucht.
    interv_search_left = -(max_disparity - offset_x_pixel) - dist_safety;
    interv_search_right = -(min_disparity - offset_x_pixel) + dist_safety;

for i = 2:-1:2
    %downsampling of images
    [img1_dwn,odd_even_xy1] = downsample(img1_rectified,i);
    [img2_dwn,odd_even_xy2] = downsample(img2_rectified,i);

    %NCC disparity map generation
    disp_map_dwn_L = sad_scanline(img1_dwn,img2_dwn,window_length,offset_x_pixel/(2^i),...
        abs(round(interv_search_right/(2^i))),abs(round(interv_search_left/(2^i))),'LR'); %1./disparity_estimation(img1_dwn,img2_dwn,...
    %round(interv_search_left/(2^i)),round(interv_search_right/(2^i)),...
    %offset_x_pixel/(2^i),window_length);
    disp_map_dwn_R = sad_scanline(img1_dwn,img2_dwn,window_length,offset_x_pixel/(2^i),...
        abs(round(interv_search_right/(2^i))),abs(round(interv_search_left/(2^i))),'RL'); %1./disparity_estimation(img1_dwn,img2_dwn,...
    %round(interv_search_left/(2^i)),round(interv_search_right/(2^i)),...
    %offset_x_pixel/(2^i),window_length);
    depth_map_dwn_L = 1./fill_disparity_map(disp_map_dwn_L,disp_map_dwn_L==0,9,5);%disp_map_dwn_L;%
    depth_map_dwn_R = 1./fill_disparity_map(disp_map_dwn_R,disp_map_dwn_R==0,9,5);%disp_map_dwn_R;%
    depth_map_up_L(:,:,i) = upsample(1/(2^i)*depth_map_dwn_L,i,odd_even_xy1);
    depth_map_up_R(:,:,i) = upsample(1/(2^i)*depth_map_dwn_R,i,odd_even_xy2);
%     figure
%     surf(depth_map_up_L(:,:,i),'LineStyle','none')
%     view(0,-90);
%     figure
%     surf(depth_map_up_R(:,:,i),'LineStyle','none')
%     view(0,-90);
    depth_map_L = depth_map_up_L(:,:,i);%min(depth_map_up,depth_map_up(:,:,i));
    depth_map_R = depth_map_up_R(:,:,i);%min(depth_map_up,depth_map_up(:,:,i));

end

% Depth Map wurde wegen der Groesse des Fensters zur SAD Berechnung an den
% Bildrändern nicht bestimmt. Darüber hinaus wurde fuer die Ermittlung der
% Depth Map nur ein Ausschnitt der rektifizierten Ansicht verwendet.
% Deshalb wird die Depth Map nun extrapoliert, indem aeussere valide Werte in
% die nicht bestimmten Bereiche kopiert werden. Aufgrund der Liniensuche
% von links nach rechts wird der rechte Rand mit der maximalen Tiefe
% versehen. Dadurch wird vermieden, dass dort aufgrund von kleinen Tiefen
% schwarze Raender entstehen.
dist_side = window_length*(2^i)/2+1;
% Maximale Tiefe bestimmen
max_depth = 1/(offset_x_pixel-interv_search_right);
% Linken Rand der Depth Map mit Werten füllen
depth_map_L(:,1:dist_side) = repmat(depth_map_L(:,dist_side+1),[1 dist_side]);
% Rechten Rand der Depth Map mit Werten füllen
% depth_map(:,end-dist_side+1:end) = repmat(depth_map(:,end-dist_side),...
%                                                             [1 dist_side]);
depth_map_L(:,end-dist_side+1:end) = ...
                               max_depth*ones(size(depth_map_L,1),dist_side);
% Oberen Rand der Depth Map mit Werten füllen
depth_map_L(1:dist_side,:) = repmat(depth_map_L(dist_side+1,:),[dist_side 1]);
% Unteren Rand der Depth Map mit Werten füllen
depth_map_L(end-dist_side+1:end,:) = repmat(depth_map_L(end-dist_side,:),...
                                                            [dist_side 1]);
% Depth map nach oben erweitern mit der bisherigen obersten Zeile
depth_map_L = [repmat(depth_map_L(1,:),d_cut_up_L,1); depth_map_L];
% Depth map nach unten erweitern mit der bisherigen untersten Zeile
depth_map_L = [depth_map_L;repmat(depth_map_L(end,:),d_cut_down_L,1)];

% Wiederhole dasselbe f�r rechtes Bild
% Linken Rand der Depth Map mit Werten fuellen
depth_map_R(:,1:dist_side) = repmat(depth_map_R(:,dist_side+1),[1 dist_side]);
% Rechten Rand der Depth Map mit Werten füllen
% depth_map(:,end-dist_side+1:end) = repmat(depth_map(:,end-dist_side),...
%                                                             [1 dist_side]);
depth_map_R(:,end-dist_side+1:end) = ...
                               max_depth*ones(size(depth_map_R,1),dist_side);
% Oberen Rand der Depth Map mit Werten füllen
depth_map_R(1:dist_side,:) = repmat(depth_map_R(dist_side+1,:),[dist_side 1]);
% Unteren Rand der Depth Map mit Werten füllen
depth_map_R(end-dist_side+1:end,:) = repmat(depth_map_R(end-dist_side,:),...
                                                            [dist_side 1]);
% Depth map nach oben erweitern mit der bisherigen obersten Zeile
depth_map_R = [repmat(depth_map_R(1,:),d_cut_up_R,1); depth_map_R];
% Depth map nach unten erweitern mit der bisherigen untersten Zeile
depth_map_R = [depth_map_R;repmat(depth_map_R(end,:),d_cut_down_R,1)];
end
% function depth_map = ...
%     depth_estimation(img1_rectified,img2_rectified,K,T,offset_x_pixel,...
%     d_cut_up,d_cut_down,min_disparity,max_disparity,varargin)
% %DEPTH_MAP 
%     window_length = 9;
%     % Minimale und maximale Disparität wurde über robuste Korrspondenzen
%     % ermittelt. dist_safety legt fest, um wie viele Pixel die tatsächliche
%     % min./max. Disparität von den ermittelten Werten abweichen darf
%     dist_safety = 20;
%     % Suchintervall festlegeben. Beispiel: Bei interv_search_left = -100
%     % und interv_search_right = 200 wird 100 Pixel nach links und 200 nach
%     % rechts gesucht.
%     interv_search_left = -(max_disparity - offset_x_pixel) - dist_safety;
%     interv_search_right = -(min_disparity - offset_x_pixel) + dist_safety;
% 
% for i = 2:-1:2
%     %downsampling of images
%     [img1_dwn,odd_even_xy1] = downsample(img1_rectified,i);
%     [img2_dwn,odd_even_xy2] = downsample(img2_rectified,i);
% 
%     %NCC disparity map generation
%     depth_map_dwn = 1./disparity_estimation(img1_dwn,img2_dwn,...
%     round(interv_search_left/(2^i)),round(interv_search_right/(2^i)),...
%     offset_x_pixel/(2^i),window_length);
%     depth_map_up(:,:,i) = upsample(1/(2^i)*depth_map_dwn,i,odd_even_xy1);
%     figure
%     surf(depth_map_up(:,:,i),'LineStyle','none')
%     view(0,-90);
%     depth_map = depth_map_up(:,:,i);%min(depth_map_up,depth_map_up(:,:,i));
% end
% 
% % Depth Map wurde wegen der Größe des Fensters zur SAD Berechnung an den
% % Bildrändern nicht bestimmt. Darüber hinaus wurde für die Ermittlung der
% % Depth Map nur ein Ausschnitt der rektifizierten Ansicht verwendet.
% % Deshalb wird die Depth Map nun extrapoliert, indem äußere valide Werte in
% % die nicht bestimmten Bereiche kopiert werden. Aufgrund der Liniensuche
% % von links nach rechts wird der rechte Rand mit der maximalen Tiefe
% % versehen. Dadurch wird vermieden, dass dort aufgrund von kleinen Tiefen
% % schwarze Ränder entstehen.
% dist_side = window_length*(2^i)/2+1;
% % Maximale Tiefe bestimmen
% max_depth = 1/(offset_x_pixel-interv_search_right);
% % Linken Rand der Depth Map mit Werten füllen
% depth_map(:,1:dist_side) = repmat(depth_map(:,dist_side+1),[1 dist_side]);
% % Rechten Rand der Depth Map mit Werten füllen
% % depth_map(:,end-dist_side+1:end) = repmat(depth_map(:,end-dist_side),...
% %                                                             [1 dist_side]);
% depth_map(:,end-dist_side+1:end) = ...
%                                max_depth*ones(size(depth_map,1),dist_side);
% % Oberen Rand der Depth Map mit Werten füllen
% depth_map(1:dist_side,:) = repmat(depth_map(dist_side+1,:),[dist_side 1]);
% % Unteren Rand der Depth Map mit Werten füllen
% depth_map(end-dist_side+1:end,:) = repmat(depth_map(end-dist_side,:),...
%                                                             [dist_side 1]);
% % Depth map nach oben erweitern mit der bisherigen obersten Zeile
% depth_map = [repmat(depth_map(1,:),d_cut_up,1); depth_map];
% % Depth map nach unten erweitern mit der bisherigen untersten Zeile
% depth_map = [depth_map;repmat(depth_map(end,:),d_cut_down,1)];
% end
