function [output_image]  = free_viewpoint(image1, image2, p)
% This function generates an image from a virtual viewpoint between two
% real images. The output image has the same size as the input images.
devMode =true;
%% Algorithmus
%% Umwandlung der Bilder in Graubilder
disp('-------------conversion to gray pictures--------------')
IGray1 = rgb_to_gray(image1);
IGray2 = rgb_to_gray(image2);

%% Intensit√§ts und Beleuchtungskorrektur
disp('-------------bias gain correction of images-----------')
IGray1_c = gain_offset_correction_cdf(IGray1);
IGray2_c = gain_offset_correction_cdf(IGray2);

%% Bilateralte Filterung
% Kanten sollen sich gut von homogenen Fl√§chen abheben k√∂nnen
%downsample
[IGray1_c, IGray2_c, odd_even] = downsample(IGray1_c,IGray2_c,1);

IGray1_bf = uint8(bfltGray(double(IGray1_c),12,100,15));
IGray2_bf = uint8(bfltGray(double(IGray2_c),12,100,15));

%upsample
[IGray1_bf] = upsample(IGray1_bf,1,odd_even);
[IGray2_bf] = upsample(IGray2_bf,1,odd_even);


%% Harris-Merkmale berechnen
disp('-------------getting feature points--------------')
% Robuste Einstellungen 'segment_length',15,'k',0.12,'min_dist',10,'N',20
Merkmale1 = harris_detektor(IGray1_bf,'segment_length',15,'k',0.15,'min_dist',20,'N',20,'do_plot',devMode);
Merkmale2 = harris_detektor(IGray2_bf,'segment_length',15,'k',0.15,'min_dist',20,'N',20,'do_plot',devMode);


%% Korrespondenzschaetzung
disp('-------------estimation correspondence points--------------')
% Robuste Einstellungen 'window_length',55,'min_corr',0.9
Korrespondenzen = punkt_korrespondenzen(IGray1_bf,IGray2_bf,Merkmale1,Merkmale2,'window_length',45,'min_corr',0.95,'do_plot',devMode);
Korrespondenzen_robust = [];
for i = 1:10
    %% Finde robuste Korrespondenzpunktpaare mit Hilfe des RANSAC-Algorithmus
    disp('-------------finding robust correspondence points--------------')
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
disp('-------------essential matrix calculation--------------')

load('calib_K.mat');

E = achtpunktalgorithmus(Korrespondenzen_robust,K);
% F = achtpunktalgorithmus(Korrespondenzen_robust);
% disp(E);
% disp(F);
[T1, R1, T2, R2, U, V] = TR_aus_E(E);

[T_cell, R_cell,T,R, d_cell, x1, x2] = ...
    rekonstruktion(T1, T2, R1, R2, Korrespondenzen_robust, K);

% load('zw_for_rec_27_08');

%% Bildrektifizierungsalgorithmus
disp('-------------rectification--------------')
[img1_rectified,img2_rectified,Tr1,Tr2,R_rect,offset_x_pixel,...
    img1_rectified_full,d_cut_up,d_cut_down] = ...
    rectification(IGray1_bf,IGray2_bf,K,T,R,'do_plot',devMode,'size_frame','valid_offset');
% figure
% imshow(img1_rectified_full)
% % save('zw_R_rect','R_rect');
% % save('comp_good')
% save('comp_28_08')
% load('comp_28_08');

% f¸r RGB-Kan‰le
for i=1:size(image1,3)
    [image1_rgb_rectified(:,:,i),image2_rgb_rectified(:,:,i),~,~,~,~,...
        image1_rgb_rectified_full(:,:,i),~,~] = ...
        rectification(image1(:,:,i),image2(:,:,i),K,T,R,'do_plot',devMode,'size_frame','valid_offset');
end

%% Disparit‰tsermittling
disp('---------disparity estimation-----------')

[min_disparity,max_disparity] = get_min_max_disparity(Korrespondenzen_robust,Tr1,Tr2);

depth_map = ...
    depth_estimation(img1_rectified,img2_rectified,K,T,offset_x_pixel,...
    d_cut_up,d_cut_down,min_disparity,max_disparity);
%% Projektion auf neues Bild
disp('---------projection-----------')
alpha = p*get_max_alpha(R);
for i=1:size(image1,3)
    [virtual_view_img(:,:,i),img_rectified_new(:,:,i)] = projection(image1(:,:,i),...
        image1_rgb_rectified_full(:,:,i),...
        depth_map,...
        K,-p,R_rect,alpha,...
        offset_x_pixel, d_cut_up,d_cut_down);
end
save('res_28_08_V3_try3')

figure
imshow(img_rectified_new);
%% Ausgabe des Free-Viewpoint Bildes
% output_image = uint8(p*image1+(1-p)*image2);
output_image = virtual_view_img;
end