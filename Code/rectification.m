function [img_L_r, img_R_r, TL, TR] = rectification (img1, img2,K,T,R,varargin)
% This function rectifies the gray images img1 and img2. The parameters,
% which must be given are EF (essential or fundamental matrix) and if EF is
% the essential matrix a camera calibratrion matrix K. Furthermore the
% correspondence_points must be given.
% Euclidean Movement is defined from Image1 to Image2

% http://www.sci.utah.edu/~gerig/CS6320-S2013/Materials/CS6320-CV-F2012-Rectification.pdf
% A compact algorithm for rectification of stereo pairs
% Andrea Fusiello1, Emanuele Trucco2, Alessandro Verri3
% http://people.scs.carleton.ca/~c_shu/Courses/comp4900d/notes/rectification.pdf

%% Input parser
p = inputParser;
valid_do_plot = @(x) islogical(x);
valid_F = @(x) size(F,1) == 3 & size(F,2) ==3;
addOptional(p,'do_plot',false,valid_do_plot);
addOptional(p,'F',zeros(0,0),valid_F);
addOptional(p,'E',zeros(0,0),valid_F);

parse(p, varargin{:});
do_plot = p.Results.do_plot;


%% Rectification Algorithum (von Paper Fusiello)
Po1 = K*[1 0 0 0; 0 1 0 0; 0 0 1 0];
Po2 = K*[R, T];
[TL, TR] = rectify_fusiello(Po1,Po2,K);

% center of left image
% Zentrieren der Bilder in Pixelkoordinaten ((0,0,1) liegt in Bildmitte)
[row_img1, col_img1] = size(img1);
[row_img2, col_img2] = size(img2);
center_row_img1 = round(row_img1/2);
center_col_img1 = round(col_img1/2);
center_row_img2 = round(row_img2/2);
center_col_img2 = round(col_img2/2);

% Berechnung & Drehung der Bildmitte, so dass Bildebenen parallel zu
% Translationsvektor liegen
center_L = [center_col_img1; center_row_img1; 1];
center_L_r = TL*center_L;
center_R = [center_col_img2; center_row_img2; 1];
center_R_r = TR*center_R;

% Vertikaler Versatz der begradigten Bilder (Y-Achse) soll nun gleich
% sein, bzw. Bilder sollen vertikal auf gleicher Höhe liegen.
dL = center_L(1:2) - center_L_r(1:2)./center_L_r(3);
dR = center_R(1:2) - center_R_r(1:2)./center_R_r(3);
dL(2) = dR(2);
[TL, TR] = rectify_fusiello(Po1,Po2,K,dL,dR);

TL_inv = inv(TL);
TR_inv = inv(TR);

% Festlegen der Größe der rektifizierten Bildausschnitte,
% Bildausschnitte sollen vorerst gleiche Größe (inkl. schwarzem Rand)
% wie Originalbilder haben
[size_L_y, size_L_x] = size(img1);
[size_R_y, size_R_x] = size(img2);

[meshX_L, meshY_L] = meshgrid(0:size_L_x-1,0:size_L_y-1);
[meshX_R, meshY_R] = meshgrid(0:size_L_x-1,0:size_L_y-1);

% Zur Interpolation werden die Bildaussschnitte zuerst in das neue
% Koordinatensystem und daraufhin zurücktransformiert. Hierfür müssen
% die Bildecken bekannt sein.
corners_L_r = TL*[0, size_L_x,size_L_x, 0;
    size_L_y, size_L_y, 0, 0;
    ones(1,4)]; % es wird der Koordinatenursprung links oben gesetzt, da TL und TR von kalibrierten Koordinaten ausgehen

corners_L_r = corners_L_r./corners_L_r(3,:);
corners_R_r = TR*[0, size_L_x,size_L_x,0;
    size_L_y, size_L_y,0,0;
    ones(1,4)];
corners_R_r = corners_R_r./corners_R_r(3,:);

% Der Bildausschnitt, der nun beide Bilder enthält, wird nun ermittelt
max_L_xy = max(ceil(corners_L_r(1:2,:)),[],2);
min_L_xy = min(floor(corners_L_r(1:2,:)),[],2);

max_R_xy = max(ceil(corners_R_r(1:2,:)),[],2);
min_R_xy = min(floor(corners_R_r(1:2,:)),[],2);

max_LR_xy = max(max_L_xy,max_R_xy);
min_LR_xy = min(min_L_xy,min_R_xy);


frame = [min_LR_xy(1), max_LR_xy(1), max_LR_xy(1), min_LR_xy(1);
    max_LR_xy(2), max_LR_xy(2), min_LR_xy(2), min_LR_xy(2);
    ones(1,4)]; %Ecke [links o., rechts o., rechts u., links.u.]

% Beide Bilder sollen die gleiche Größe haben,
% Um interpolieren zu können wird ein meshgrid benötigt, das
% zurücktransformiert werden kann
[meshX_L_r,meshY_L_r] = meshgrid(min_LR_xy(1):max_LR_xy(1)-1,...
    min_LR_xy(2):max_LR_xy(2)-1);
meshZ_L_r = ones(size(meshX_L_r));
[meshX_R_r,meshY_R_r] = meshgrid(min_LR_xy(1):max_LR_xy(1)-1,...
    min_LR_xy(2):max_LR_xy(2)-1);
meshZ_R_r = ones(size(meshX_R_r));

% Rücktransformation der Meshgrids zur Interpolation
meshZ_L_dr = TL_inv(3,1)*meshX_L_r + TL_inv(3,2)*meshY_L_r ...
    + TL_inv(3,3)*meshZ_L_r;
meshX_L_dr = (TL_inv(1,1)*meshX_L_r + TL_inv(1,2)*meshY_L_r...
    + TL_inv(1,3)*meshZ_L_r)./meshZ_L_dr;
meshY_L_dr = (TL_inv(2,1)*meshX_L_r + TL_inv(2,2)*meshY_L_r...
    + TL_inv(2,3)*meshZ_L_r)./meshZ_L_dr;

meshZ_R_dr = TR_inv(3,1)*meshX_R_r + TR_inv(3,2)*meshY_R_r...
    + TR_inv(3,3)*meshZ_R_r;
meshX_R_dr = (TR_inv(1,1)*meshX_R_r + TR_inv(1,2)*meshY_R_r...
    + TR_inv(1,3)*meshZ_R_r)./meshZ_R_dr;
meshY_R_dr = (TR_inv(2,1)*meshX_R_r + TR_inv(2,2)*meshY_R_r...
    + TR_inv(2,3)*meshZ_R_r)./meshZ_R_dr;

img_L_r = uint8(interp2(meshX_L,meshY_L,double(img1),...
    meshX_L_dr,meshY_L_dr,'linear',NaN));
img_R_r = uint8(interp2(meshX_R,meshY_R,double(img2),...
    meshX_R_dr,meshY_R_dr,'linear',NaN));

if(do_plot)
    figure
    subplot(1,2,1)
    imshow(uint8(img_L_r));
    subplot(1,2,2)
    imshow(uint8(img_R_r));
    
    figure
    im1 = imshow(uint8(img_L_r));hold all;
    im2 = imshow(uint8(img_R_r));
    im1.AlphaData = 0.5;
    im2.AlphaData = 0.5;
end

%% TO BE REMOVED

% version = 12;
% % Rektifizierung anhand Epipolen (x-Achse soll auf jeweilige Epipole gelegt
% % werden
% if(version == 1)
%     e1 = -R'*T;
%     e1 = e1./e1(end);
%     e2 = -T./T(end); % x-Achse in Bild 2 soll in gleiche Richtung wie x-Achse in Bild 1 zeigen
%
%
%     r11 = e1./sqrt(sum(e1.^2));
%     r21 = [-e1(2); e1(1); 0]./sqrt(e1(1)^2+e1(2)^2);
%     r31 = cross(r11,r21);
%
%     R_rect1 = [r11, r21, r31];
%
%     r12 = e2./sqrt(sum(e2.^2));
%     r22 = [-e2(2); e2(1); 0]./sqrt(e2(1)^2+e2(2)^2);
%     r32 = cross(r12,r22);
%
%     R_rect2 = [r12'; r22'; r32'];%R_rect1*R';
%
%     T1 = K*R_rect1*inv(K);
%     T2 = K*R_rect1*R'*inv(K);
%
% elseif(version ==2)
%     [U,SIGMA,V] = svd(F);
%     f1 = V(:,3)./V(3,3);
%     %     f1(2) = 0;
%     f2 = -U(:,3)./U(3,3);
%     %     f2(2) = 0;
%
%     r11 = f1./sqrt(sum(f1.^2));
%     r21 = [-f1(2); f1(1); 0]./sqrt(f1(1)^2+f1(2)^2);
%     r31 = cross(r11,r21);
%
%     R_rect1 = [r11, r21, r31];
%
%     r12 = f2./sqrt(sum(f2.^2));
%     r22 = [-f2(2); f2(1); 0]./sqrt(f2(1)^2+f2(2)^2);
%     r32 = cross(r12,r22);
%     R_rect2 = [r12, r22, r32];%R_rect1*R';
%
%     T1 = R_rect1;
%     T2 = R_rect2';
% end

% T1_inv = inv(T1);
% T2_inv = inv(T2);
%
% % % Zentrieren der Bilder in Pixelkoordinaten ((0,0,1) liegt in Bildmitte)
% % [row_img1, col_img1] = size(img1);
% % [row_img2, col_img2] = size(img2);
% % center_row_img1 = round(row_img1/2);
% % center_col_img1 = round(col_img1/2);
% % center_row_img2 = round(row_img2/2);
% % center_col_img2 = round(col_img2/2);
%
% [X_img1,Y_img1] = meshgrid(-center_col_img1:center_col_img1-1,...
%     -center_row_img1:center_row_img1-1);
% [X_img2,Y_img2] = meshgrid((-center_col_img2:center_col_img2-1),...
%     (-center_row_img2:center_row_img2-1));
%
% % Eckpunkte vom img1 und img2 berechnen & ins neue Koordinatensystem
% % transformieren
% corners_img1 = [-center_col_img1, center_col_img1,center_col_img1,-center_col_img1;
%     center_row_img1, center_row_img1, -center_row_img1, -center_row_img1;
%     ones(1,4)];
% corners_img2 = [-center_col_img2, center_col_img2,center_col_img2,-center_col_img2;
%     center_row_img2, center_row_img2, -center_row_img2, -center_row_img2;
%     ones(1,4)];
% % Plotting der Corner-Frames
% figure
% plot3([corners_img1(1,:) corners_img1(1,1)],[corners_img1(2,:) corners_img1(2,1)],[corners_img1(3,:) corners_img1(3,1)],'b');
% hold all
% grid
% corners2 = K*R'*(inv(K)*corners_img1-T);
% plot3([corners2(1,:) corners2(1,1)],[corners2(2,:) corners2(2,1)],[corners2(3,:) corners2(3,1)],'r');
% o1 = zeros(3,1);
% o2 = -K*R'*T;
% plot3(o1(1),o1(2),o1(3),'xb')
% plot3(o2(1),o2(2),o2(3),'xr')
% %plot3([0 f1(1)],[0 f1(2)],[0 f1(3)]);
% % plot3([0 f2(1)],[0 f2(2)],[0 f2(3)]);
%
% corners_img1_r = T1*corners_img1;
% corners_img2_r = T2*corners_img2-K*R'*T;
% plot3([corners_img1_r(1,:) corners_img1_r(1,1)],[corners_img1_r(2,:) corners_img1_r(2,1)],[corners_img1_r(3,:) corners_img1_r(3,1)],'-.b');
% plot3([corners_img2_r(1,:) corners_img2_r(1,1)],[corners_img2_r(2,:) corners_img2_r(2,1)],[corners_img2_r(3,:) corners_img2_r(3,1)],'-.r');
% % axis equal
% xlabel('X'); ylabel('Y'); zlabel('Z');
% % rektifizierte Bildecken (gerundet)
% corners_img1_r = T1*corners_img1;
% corners_img2_r = T2*corners_img2;
% corners_img1_r = round(corners_img1_r./corners_img1_r(3,:));
% corners_img2_r = round(corners_img2_r./corners_img2_r(3,:));
%
%
% % erzeugen eines Meshgrids für das rektifizierte Bild
% [X_img1_r, Y_img1_r] = meshgrid(min(corners_img1_r(1,:)):max(corners_img1_r(1,:)),...
%     min(corners_img1_r(2,:)):max(corners_img1_r(2,:)));
% [X_img2_r, Y_img2_r] = meshgrid((min(corners_img2_r(1,:)):max(corners_img2_r(1,:))),...
%     min(corners_img2_r(2,:)):max(corners_img2_r(2,:)));
% Z_img1_r = ones(size(X_img1_r));
% Z_img2_r = ones(size(X_img2_r));
%
% % Rücktransformation der Meshgrids zur Interpolation
% Z_img1_dr = T1_inv(3,1)*X_img1_r + T1_inv(3,2)*Y_img1_r + T1_inv(3,3)*Z_img1_r;
% X_img1_dr = (T1_inv(1,1)*X_img1_r + T1_inv(1,2)*Y_img1_r + T1_inv(1,3)*Z_img1_r)./Z_img1_dr;
% Y_img1_dr = (T1_inv(2,1)*X_img1_r + T1_inv(2,2)*Y_img1_r + T1_inv(2,3)*Z_img1_r)./Z_img1_dr;
%
% Z_img2_dr = T2_inv(3,1)*X_img2_r + T2_inv(3,2)*Y_img2_r + T2_inv(3,3)*Z_img2_r;
% X_img2_dr = (T2_inv(1,1)*X_img2_r + T2_inv(1,2)*Y_img2_r + T2_inv(1,3)*Z_img2_r)./Z_img2_dr;
% Y_img2_dr = (T2_inv(2,1)*X_img2_r + T2_inv(2,2)*Y_img2_r + T2_inv(2,3)*Z_img2_r)./Z_img2_dr;
%
% min_X_img1_r = min(corners_img1_r(1,:));
% min_X_img2_r = min(corners_img2_r(1,:));
% min_Y_img1_r = min(corners_img1_r(2,:));
% min_Y_img2_r = min(corners_img2_r(2,:));
% max_X_img1_r = max(corners_img1_r(1,:));
% max_X_img2_r = max(corners_img2_r(1,:));
% max_Y_img1_r = max(corners_img1_r(2,:));
% max_Y_img2_r = max(corners_img2_r(2,:));
%
%
% % Zwei einzelne Bilder
% img_1u_r = zeros(length(min_Y_img1_r:max_Y_img1_r),length(min_X_img1_r:max_X_img1_r));
% img_1u_r((min_Y_img1_r:max_Y_img1_r)+1-min_Y_img1_r,...
%     (min_X_img1_r:max_X_img1_r)+1-min_X_img1_r) ...
%     = interp2(X_img1,Y_img1,double(img1),X_img1_dr,Y_img1_dr,'linear',NaN);
% img_2u_r = zeros(length(min_Y_img2_r:max_Y_img2_r),length(min_X_img2_r:max_X_img2_r));
% img_2u_r((min_Y_img2_r:max_Y_img2_r)+1-min_Y_img2_r,...
%     (min_X_img2_r:max_X_img2_r)+1-min_X_img2_r) ...
%     = interp2(X_img2,Y_img2,double(img2),X_img2_dr,Y_img2_dr,'linear',NaN);
%
% % Einzelne Bilder übereinander legen (ausgehend vom Zentrum der Bilder)
% [size_y_im1, size_x_im1] = size(img_1u_r);
% [size_y_im2, size_x_im2] = size(img_2u_r);
%
% center_y_frame = max(ceil(size_y_im1/2),ceil(size_y_im2/2));
% center_x_frame = max(ceil(size_x_im1/2),ceil(size_x_im2/2));
% max_dist_y_frame1 = floor(size_y_im1/2);
% if(mod(size_y_im1,2))
%     % nicht durch 2 teilbar
%     min_dist_y_frame1 = -floor(size_y_im1/2);
% else
%     % durch 2 teilbar
%     min_dist_y_frame1 = -floor(size_y_im1/2)+1;
% end
% max_dist_x_frame1 = floor(size_x_im1/2);
% if(mod(size_x_im1,2))
%     min_dist_x_frame1 = -floor(size_x_im1/2);
% else
%     min_dist_x_frame1 = -floor(size_x_im1/2)+1;
% end
%
% max_dist_y_frame2 = floor(size_y_im2/2);
% if(mod(size_y_im2,2))
%     min_dist_y_frame2 = -floor(size_y_im2/2);
% else
%     min_dist_y_frame2 = -floor(size_y_im2/2)+1;
% end
% max_dist_x_frame2 = floor(size_x_im2/2);
% if(mod(size_x_im2,2))
%     min_dist_x_frame2 = -floor(size_x_im2/2);
% else
%     min_dist_x_frame2 = -floor(size_x_im2/2)+1;
% end
%
% frame1 = zeros(max(size_y_im1,size_y_im2),max(size_y_im2,size_x_im2));
% frame2 = frame1;
%
% frame1((min_dist_y_frame1:max_dist_y_frame1)+center_y_frame,...
%     (min_dist_x_frame1:max_dist_x_frame1)+center_x_frame) = img_1u_r;
% frame2((min_dist_y_frame2:max_dist_y_frame2)+center_y_frame,...
%     (min_dist_x_frame2:max_dist_x_frame2)+center_x_frame) = img_2u_r;
%
% figure
% subplot(1,2,1)
% imshow(uint8(img_1u_r));
% subplot(1,2,2)
% imshow(uint8(img_2u_r));
%
% figure
% im1 = imshow(uint8(frame1));hold all;
% im2 = imshow(uint8(frame2));
% im1.AlphaData = 0.5;
% im2.AlphaData = 0.5;

end
