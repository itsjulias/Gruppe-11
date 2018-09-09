function depth_map = ...
    depth_estimation(img1_rectified,img2_rectified,offset_x_pixel,...
    d_cut_down,d_cut_up,min_disparity,max_disparity,log_NaN_L,log_NaN_R,mean_disparity,L_R_base,varargin)
% DEPTH_ESTIMATION Die Tiefen werden mit Hilfe den downgesampelten Bilder
% berechnet. Zunächst wird die geringste Auflösung (Layer 3, 2x
% Downsampling) verwendet. Mit Hilfe von SAD werden zunächst grob die
% Disparitäten bestimmt. Gleichzeitig wird die Varianz jedes
% Fensterausschnitts verwendet, um homogene Flächen zu bestimmen. Die
% Disparitäten der homogenen Flächen werden nicht genauer bestimmt, da dies
% auch nur schwer möglich ist. Die restlichen Stellen werden hingegen in
% der nächst höheren Auflösung (Layer 2) genauer definiert. Um Rechenzeit zu sparen
% werden die grob bestimmten Disparitäten der niedrigeren Auflösung
% verwendet, um für jedes Pixel ein (kleineres) Liniensuchintervall zu
% erhalten (=> Erneute Suche mit refine_disp.m).
% Für kleine Disparitäten kann auch ein gröbere Depth Map ausreichen und
% Rechenzeit gespart werden. (Fall mean_disparity <= 1500)

% Quelle:
% Fast Computation of Dense and Reliable Depth Maps from Stereo Images
% M. Tornow, M. Grasshoff, N. Nguyen, A. Al-Hamadi and B. Michaelis

    % Minimale und maximale Disparitüt wurde Ã¼ber robuste Korrspondenzen
    % ermittelt. dist_safety legt fest, um wie viele Pixel die tatsÃ¤chliche
    % min./max. DisparitÃ¤t von den ermittelten Werten abweichen darf
    dist_safety = 20;
    % Suchintervall festlegeben. Beispiel: Bei interv_search_left = -100
    % und interv_search_right = 200 wird 100 Pixel nach links und 200 nach
    % rechts gesucht.
    interv_search_left = -(max_disparity - offset_x_pixel) - dist_safety;
    interv_search_right = -(min_disparity - offset_x_pixel) + dist_safety;        
    
% Bei großen Disparitäten: Bestimme Disparitäten auch in der nächst höheren
% Auflösung genauer.
if mean_disparity > 1500 % true für Bilder L1/R1
    %% Beginne mit dem 3. Layer (d.h. zweimal downsampling)
    disp('(Using Refinement)')
    i = 2;
    window_length(1) = 5;
    [img2_dwn, odd_even_xy2] = downsample(img2_rectified,i,'no');
    [img1_dwn, odd_even_xy] = downsample(img1_rectified,i,'no');
    
    if strcmp(L_R_base,'R_base')
        % Falls rechte Ansicht die Ausgangsansicht ist.
        img = img1_dwn;
        img1_dwn = img2_dwn;
        img2_dwn = img;
        interv = interv_search_left;
        interv_search_left = -interv_search_right;
        interv_search_right = -interv;
        odd_even_xy = odd_even_xy2;
    end
    
    % Disparität wird mit den downgesampelten Bildern berechnet.
    [disp_map_3,hom_map] = disparity_estimation(img1_dwn,img2_dwn,...
    round(interv_search_left/(2^i)),round(interv_search_right/(2^i)),...
    offset_x_pixel/(2^i),window_length(1),L_R_base);

%     disp_surf(:,:,1) = disp_map_3;
%     figure
%     surf(disp_surf(:,:,1),'LineStyle','none')
%     view(0,-90);
    
    %% Median Filter 
    % Homogene Flächen sollen möglichst keine Ausreißer enthalten.
    disp_map_filtered_3 = medfilt2(disp_map_3,[9 9]);

%     disp_surf(:,:,1) = disp_map_filtered_3;
%     figure
%     surf(disp_surf(:,:,1),'LineStyle','none')
%     view(0,-90);
    
    disp_map_filtered_3_up = upsample((2^i)*disp_map_filtered_3,1,odd_even_xy(2,:));
    
    %% Layer 2: Disparitätsermittlung  
    disp_map_2 = upsample((2^i)*disp_map_3,1,odd_even_xy(2,:));
    hom_map = upsample(hom_map,1,odd_even_xy(2,:));
    
    % Ersetze mögliche NaNs durch upsamling durch 0
    hom_map(isnan(hom_map)) = 0;
    
    % Min./Max. Disparität für jedes Pixel der Disparity Map bestimmen.
    win_len_minmax = 9;
    min_max_disp = [];
    min_max_disp(:,:,1) = movmin(movmin(disp_map_2,win_len_minmax),win_len_minmax,2);
    min_max_disp(:,:,2) = movmax(movmax(disp_map_2,win_len_minmax),win_len_minmax,2);
    
%     disp_surf = [];
%     disp_surf(:,:,1) = min_max_disp(:,:,2)-min_max_disp(:,:,1);
%     figure
%     surf(disp_surf(:,:,1),'LineStyle','none')
%     view(0,-90);
    
    % Neues Suchintervall für jedes Pixel festlegen
    % (vgl. Rechnung zu Beginn dieser Function)
    % Individueller Suchbereich zu jedem Pixel muss innerhalb des generell
    % definierten maximalen Suchbereichs liegen, der durch
    % interv_search_left und interv_search_right definiert ist.
    interv_search_left_matrix = -min(round(min_max_disp(:,:,2) - offset_x_pixel),-interv_search_left);
    interv_search_right_matrix = min(-round(min_max_disp(:,:,1) - offset_x_pixel),interv_search_right);

    
    %% Bestimme jeweils Layer 2
    i = 1;
    [img2_dwn, odd_even_xy2] = downsample(img2_rectified,i,'no');
    [img1_dwn, odd_even_xy] = downsample(img1_rectified,i,'no');
    
    if strcmp(L_R_base,'R_base')
        % Falls rechte Ansicht die Ausgangsansicht ist.
        img = img1_dwn;
        img1_dwn = img2_dwn;
        img2_dwn = img;
        interv = interv_search_left_matrix;
        interv_search_left_matrix = -interv_search_right_matrix;
        interv_search_right_matrix = -interv;
        odd_even_xy = odd_even_xy2;
    end

    %% Refine
    window_length(2) = 17;
    % Bestimme anhand des für jedes Pixel individuell festgelegten
    % Liniensuchbereichs erneut die Disparitäten.
    disp_map_refine = refine_disp(img1_dwn,img2_dwn,hom_map,...
        round(interv_search_left_matrix./(2^i)),...
        round(interv_search_right_matrix./(2^i)),...
            offset_x_pixel/(2^i),window_length(2),L_R_base);
    
%     disp_surf = [];
%     disp_surf(:,:,1) = (2^i)*disp_map_refine;
%     figure
%     surf(disp_surf(:,:,1),'LineStyle','none')
%     view(0,-90);


    %% Merge
    % Bestimme Stellen, an denen disp_map_refine nicht definiert ist und
    % fülle diese Stellen mit der gefilterten Depth Map
    disp_merge = (2^i)*disp_map_refine;
    disp_merge(isnan(disp_map_refine)) = disp_map_filtered_3_up(isnan(disp_map_refine));
    
%     disp_surf = [];
%     disp_surf(:,:,1) = disp_merge;
%     figure
%     surf(disp_surf(:,:,1),'LineStyle','none')
%     view(0,-90);
        
    disp_map_filtered = continuity_filter(disp_merge,3,9); %,3
    
%     disp_surf = [];
%     disp_surf(:,:,1) = disp_map_filtered;
%     figure
%     surf(disp_surf(:,:,1),'LineStyle','none')
%     view(0,-90);
    
    % Depth map generation
    depth_map_dwn = 1./disp_map_filtered;
    depth_map_up(:,:,1) = upsample(depth_map_dwn,i,odd_even_xy);
    
else % true für Bilder L2/R2
    %% Nur 3. Layer berechnen (d.h. zweimal downsampling)
    disp('(Only using Layer 3)')
    i = 2;
    window_length(1) = 11;
    [img2_dwn, odd_even_xy2] = downsample(img2_rectified,i,'no');
    [img1_dwn, odd_even_xy] = downsample(img1_rectified,i,'no');
    
    if strcmp(L_R_base,'R_base')
        % Falls rechte Ansicht die Ausgangsansicht ist.
        img = img1_dwn;
        img1_dwn = img2_dwn;
        img2_dwn = img;
        interv = interv_search_left;
        interv_search_left = -interv_search_right;
        interv_search_right = -interv;
        odd_even_xy = odd_even_xy2;
    end
    
    % Disparität wird mit den downgesampelten Bildern berechnet.
    [disp_map_3,~] = disparity_estimation(img1_dwn,img2_dwn,...
    round(interv_search_left/(2^i)),round(interv_search_right/(2^i)),...
    offset_x_pixel/(2^i),window_length(1),L_R_base);
    % Layer 3 reicht aus um ausreichende Depth Map zu erhalten.
    disp_map_filtered = continuity_filter(disp_map_3,3,7);
    disp_map_up = upsample((2^i)*disp_map_filtered,2,odd_even_xy);
    
%     disp_surf = [];
%     disp_surf(:,:,1) = disp_map_up;
%     figure
%     surf(disp_surf(:,:,1),'LineStyle','none')
%     view(0,-90);

    depth_map_up(:,:,1) = 1./disp_map_up;

%     figure
%     surf(depth_map_up(:,:,1),'LineStyle','none')
%     view(0,-90);

    window_length(2) = window_length(1);
end

    depth_map = depth_map_up(:,:,1);

%% Expand
% Depth Map wurde wegen der Größe des Fensters zur SAD Berechnung an den
% Bildrändern nicht bestimmt. Darüber hinaus wurde für die Ermittlung der
% Depth Map nur ein Ausschnitt der rektifizierten Ansicht verwendet.
% Deshalb wird die Depth Map nun extrapoliert, indem äußere valide Werte in
% die nicht bestimmten Bereiche kopiert werden. Aufgrund der Liniensuche
% von links nach rechts wird der rechte Rand mit der maximalen Tiefe
% versehen. Dadurch wird vermieden, dass dort aufgrund von kleinen Tiefen
% schwarze Ränder entstehen.
i = 2;
dist_side = window_length(2)*(2^i)/2+1;
% Maximale Tiefe bestimmen
max_depth = 1/(offset_x_pixel-interv_search_right);
% Linken Rand der Depth Map mit Werten füllen
depth_map(:,1:dist_side) = repmat(depth_map(:,dist_side+1),[1 dist_side]);
% Rechten Rand der Depth Map mit Werten füllen
depth_map(:,end-dist_side+1:end) = ...
                               max_depth*ones(size(depth_map,1),dist_side);
% Oberen Rand der Depth Map mit Werten füllen
depth_map(1:dist_side,:) = repmat(depth_map(dist_side+1,:),[dist_side 1]);
% Unteren Rand der Depth Map mit Werten füllen
depth_map(end-dist_side+1:end,:) = repmat(depth_map(end-dist_side,:),...
                                                            [dist_side 1]);
% Depth map nach oben erweitern mit der bisherigen obersten Zeile
depth_map = [repmat(depth_map(1,:),d_cut_down,1); depth_map];
% Depth map nach unten erweitern mit der bisherigen untersten Zeile
depth_map = [depth_map;repmat(depth_map(end,:),d_cut_up,1)];

% Rektifiziertes Bild ist nicht rechteckig => Übernehme die nicht gültigen
% Bereiche der Rektifizierung
    if strcmp(L_R_base,'R_base')
        depth_map(log_NaN_R) = NaN;
    else
        depth_map(log_NaN_L) = NaN;
    end

end
