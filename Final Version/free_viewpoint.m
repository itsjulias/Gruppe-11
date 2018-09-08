function [output_image]  = free_viewpoint(image1, image2, p)
% This function generates an image from a virtual viewpoint between two
% real images. The output image has the same size as the input images.
devMode =false;
if(p == 1)
    output_image = image2;
elseif(p==0)
    output_image = image1;
else
    %% Algorithmus
    %% Umwandlung der Bilder in Graubilder
    disp('-------------CONVERSION TO GRAY IMAGES--------------')
    IGray1 = rgb_to_gray(image1);
    IGray2 = rgb_to_gray(image2);
    
    %% IntensitÃ¤ts und Beleuchtungskorrektur
    disp('-------------BIAS GAIN OFFSET CORRECTION-----------')
    IGray1_c = gain_offset_correction_cdf(IGray1);
    IGray2_c = gain_offset_correction_cdf(IGray2);
    
    %% Bilateralte Filterung
    % Kanten sollen sich gut von homogenen FlÃ¤chen abheben kÃ¶nnen
    %downsample
    [IGray1_c, odd_even1] = downsample(IGray1_c,1);
    [IGray2_c, odd_even2] = downsample(IGray2_c,1);
    
    IGray1_bf = uint8(bfltGray(double(IGray1_c),12,100,15));
    IGray2_bf = uint8(bfltGray(double(IGray2_c),12,100,15));
    
    %upsample
    [IGray1_bf] = upsample(IGray1_bf,1,odd_even1);
    [IGray2_bf] = upsample(IGray2_bf,1,odd_even2);
    
    
    %% Harris-Merkmale berechnen
    disp('-------------GETTING FEATURE POINTS--------------')
    Merkmale1 = harris_detektor(IGray1_bf,'segment_length',15,'k',0.15,'min_dist',20,'N',20,'do_plot',devMode);
    Merkmale2 = harris_detektor(IGray2_bf,'segment_length',15,'k',0.15,'min_dist',20,'N',20,'do_plot',devMode);
    
    
    %% Korrespondenzschaetzung
    disp('-------------CORRESPONDENCE POINT ESTIAMATION--------------')
    % Robuste Einstellungen 'window_length',55,'min_corr',0.9
    Korrespondenzen = punkt_korrespondenzen(IGray1_bf,IGray2_bf,Merkmale1,Merkmale2,'window_length',45,'min_corr',0.95,'do_plot',devMode);
    Korrespondenzen_robust = [];

    %% Finde robuste Korrespondenzpunktpaare mit Hilfe des RANSAC-Algorithmus
    disp('-------------FINDING ROBUST CORRESPONDENCE POINTS--------------')
    for i = 1:10
        disp(sprintf('-> Iteration %i of 10',i) );
        Korrespondenzen_robust = [Korrespondenzen_robust,...
            F_ransac(Korrespondenzen,'epsilon',0.7,'p',0.999,'tolerance', 0.01)];
        
        % Zeige die robusten Korrespondenzpunktpaare
        if(devMode)
            disp('-------------plot robust correspondence points--------------')
            figure
            im1 = imshow(IGray1_bf);
            hold all;
            im2 = imshow(IGray2_bf);
            im1.AlphaData = 0.5;
            im2.AlphaData = 0.5;
            plot(Korrespondenzen_robust(1,:),Korrespondenzen_robust(2,:),'*r');
            plot(Korrespondenzen_robust(3,:),Korrespondenzen_robust(4,:),'*g');
            for i = 1:length(Korrespondenzen_robust(1,:))
                plot([Korrespondenzen_robust(1,i) Korrespondenzen_robust(3,i)],[Korrespondenzen_robust(2,i) Korrespondenzen_robust(4,i)],'b');
            end
        end
    end
    Korrespondenzen_robust = unique(Korrespondenzen_robust','rows')';
    
    %% Berechne die Essentielle Matrix
    % Kamerakalibrierungsmatrix ist in KK.mat enthalten und wurde mit in der
    % Video-Vorlesung angegebenen Toolbox bestimmt.
    disp('-------------ESSENTIAL MATRIX ESTIMATION--------------')
    
    load('calib_K.mat');
    E = achtpunktalgorithmus(Korrespondenzen_robust,K);
    [T1, R1, T2, R2, U, V] = TR_aus_E(E);
    [T_cell, R_cell,T,R, d_cell, x1, x2] = ...
        rekonstruktion(T1, T2, R1, R2, Korrespondenzen_robust, K);
    
    %% Bildrektifizierungsalgorithmus
    disp('-------------RECTIFICATION--------------')
    
    [img1_rectified,img2_rectified,Tr1,Tr2,R_rect,offset_x_pixel,...
    img1_rectified_full,img2_rectified_full,d_cut_down_L,d_cut_up_L,...
    d_cut_down_R,d_cut_up_R,min_corners_L_y,log_NaN_L,log_NaN_R] = ...
    rectification(IGray1,IGray2,K,T,R);
    
    % fuer RGB-Kanaele
    for i=1:size(image1,3)
        [image1_rgb_rectified(:,:,i),image2_rgb_rectified(:,:,i),~,~,~,~,...
            image1_rgb_rectified_full(:,:,i),image2_rgb_rectified_full(:,:,i),~,~,~,~,~,~,~] = ...
            rectification(image1(:,:,i),image2(:,:,i),K,T,R);
    end
    

    %% Disparitaetsermittlung
    disp('---------DEPTH ESTIMATION-----------')
    
    [min_disparity,max_disparity,mean_disparity] = get_stat_disparity(Korrespondenzen_robust,Tr1,Tr2);
    
    disp('-> Calculating left depth map')
    depth_map_L = ...
        depth_estimation(img1_rectified,img2_rectified,offset_x_pixel,...
        d_cut_down_L,d_cut_up_L,min_disparity,max_disparity,...
        log_NaN_L,log_NaN_R,mean_disparity,'L_base');
    
    disp('-> Calculating right depth map')
    depth_map_R = ...
        depth_estimation(img1_rectified,img2_rectified,offset_x_pixel,...
        d_cut_down_R,d_cut_up_R,min_disparity,max_disparity,log_NaN_L,...
        log_NaN_R,mean_disparity,'R_base');
    

    %% Projektion auf neues Bild
    disp('---------PROJECTION-----------')
    disp('-> Projection of R-Channel')
    [virtual_view_img(:,:,1),img_rectified_new(:,:,1),ScatInterp] = projection(IGray1,...
    image1_rgb_rectified_full(:,:,1),image2_rgb_rectified_full(:,:,1),...
    {depth_map_L,depth_map_R},K,R,p,R_rect,...
    offset_x_pixel,min_corners_L_y);
    
    for i = 2:size(image1,3)
        if i == 2
            disp('-> Projection of G-Channel')
        else
            disp('-> Projection of B-Channel')
        end
        [virtual_view_img(:,:,i),img_rectified_new(:,:,i),~] = projection(IGray1,...
            image1_rgb_rectified_full(:,:,i),image2_rgb_rectified_full(:,:,i),...
            {},K,R,p,R_rect,...
            offset_x_pixel,min_corners_L_y,ScatInterp);
    end
    
    
    %% Ausgabe des Free-Viewpoint Bildes
    output_image = virtual_view_img;
end

end
