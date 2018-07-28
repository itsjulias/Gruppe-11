function [img_rect1, img_rect2, R_rect] = rectification (img1, img2, T, K, R,correspondence_points, EF)
% This function rectifies the gray images img1 and img2. The parameters,
% which must be given are EF (essential or fundamental matrix) and if EF is
% the essential matrix a camera calibratrion matrix K. Furthermore the
% correspondence_points must be given.
% Euclidean Movement is defined from Image1 to Image2

% http://www.sci.utah.edu/~gerig/CS6320-S2013/Materials/CS6320-CV-F2012-Rectification.pdf
% A compact algorithm for rectification of stereo pairs
% Andrea Fusiello1, Emanuele Trucco2, Alessandro Verri3
% http://people.scs.carleton.ca/~c_shu/Courses/comp4900d/notes/rectification.pdf


% Rotationsmatrix zur Bildrefktifizierung (auf Basis des
% Translationsvektors T -> zeigt von Cam2 (rechts) nach Cam1 (linkts))
T = -R'*T; % Cam1 (links) ist in Weltkoordinaten (0,0,0)
r1 = T/sqrt(sum(T.^2));
r2 = [-T(2); T(1); 0]./sqrt(T(1)^2+T(2)^2);
r3 = cross(r1,r2);

R_rect = [r1'; r2'; r3'];

% Kamerazentrum Cam2 in Pixelkoordinatnen
Cam2_pix = round(K*R_rect*T);

% neue Projektionsmatrix
pml1 = K*[R_rect, zeros(3,1)];
pmr1 = K*[R_rect, -R_rect*T];

%Transformationsvorschrift zur Bildrektifizierung
T1 = K*R_rect*inv(K);
T2 = K*R_rect*inv(K*R');%K*R_rect*(R')*inv(K);

T1_inv = inv(T1);
T2_inv = inv(T2);

[row_img1, col_img1] = size(img1);
[row_img2, col_img2] = size(img2);
% Zentrieren der Bilder in Pixelkoordinaten ((0,0,1) liegt in Bildmitte)
% tbd: rounding wenn Bild 3001x2001 pixel hat
center_row_img1 = round(row_img1/2);
center_col_img1 = round(col_img1/2);
center_row_img2 = round(row_img2/2);
center_col_img2 = round(col_img2/2);

[X_img1,Y_img1] = meshgrid(-center_col_img1:center_col_img1-1,...
                            -center_row_img1:center_row_img1-1);
[X_img2,Y_img2] = meshgrid((-center_col_img2:center_col_img2-1)+Cam2_pix(1),...
                           (-center_row_img2:center_row_img2-1)+Cam2_pix(2));

% Eckpunkte vom img1 und img2 berechnen & ins neue Koordinatensystem
% transformieren
corners_img1 = [-center_col_img1, center_col_img1,-center_col_img1, center_col_img1;
                center_row_img1, center_row_img1, -center_row_img1, -center_row_img1;
                ones(1,4)];
corners_img2 = [-center_col_img2, center_col_img2,-center_col_img2, center_col_img2;
                center_row_img2, center_row_img2, -center_row_img2, -center_row_img2;
                ones(1,4)];        
% rektifizierte Bildecken (gerundet)
corners_img1_r = T1*corners_img1;
corners_img2_r = T2*corners_img2;
corners_img1_r = round(corners_img1_r./corners_img1_r(3,:));
corners_img2_r = round(corners_img2_r./corners_img2_r(3,:));

% erzeugen eines Meshgrids für das rektifizierte Bild
[X_img1_r, Y_img1_r] = meshgrid(min(corners_img1_r(1,:)):max(corners_img1_r(1,:)),...
                                min(corners_img1_r(2,:)):max(corners_img1_r(2,:)));
[X_img2_r, Y_img2_r] = meshgrid((min(corners_img2_r(1,:)):max(corners_img2_r(1,:)))+Cam2_pix(1),...
                            min(corners_img2_r(2,:)):max(corners_img2_r(2,:))+Cam2_pix(2));
Z_img1_r = ones(size(X_img1_r));
Z_img2_r = ones(size(X_img2_r));

% Rücktransformation der Meshgrids zur Interpolation
Z_img1_dr = T1_inv(3,1)*X_img1_r + T1_inv(3,2)*Y_img1_r + T1_inv(3,3)*Z_img1_r;
X_img1_dr = (T1_inv(1,1)*X_img1_r + T1_inv(1,2)*Y_img1_r + T1_inv(1,3)*Z_img1_r)./Z_img1_dr;
Y_img1_dr = (T1_inv(2,1)*X_img1_r + T1_inv(2,2)*Y_img1_r + T1_inv(2,3)*Z_img1_r)./Z_img1_dr;

Z_img2_dr = T2_inv(3,1)*X_img2_r + T2_inv(3,2)*Y_img2_r + T2_inv(3,3)*Z_img2_r;
X_img2_dr = (T2_inv(1,1)*X_img2_r + T2_inv(1,2)*Y_img2_r + T2_inv(1,3)*Z_img2_r)./Z_img2_dr;
Y_img2_dr = (T2_inv(2,1)*X_img2_r + T2_inv(2,2)*Y_img2_r + T2_inv(2,3)*Z_img2_r)./Z_img2_dr;

min_X_img1_r = min(corners_img1_r(1,:));
min_X_img2_r = min(corners_img2_r(1,:));
min_Y_img1_r = min(corners_img1_r(2,:));
min_Y_img2_r = min(corners_img2_r(2,:));
max_X_img1_r = max(corners_img1_r(1,:));
max_X_img2_r = max(corners_img2_r(1,:));
max_Y_img1_r = max(corners_img1_r(2,:));
max_Y_img2_r = max(corners_img2_r(2,:));

min_Y_img_r = min(min_Y_img1_r,min_Y_img2_r);
max_Y_img_r = max(max_Y_img1_r,max_Y_img2_r);
min_X_img_r = min(min_X_img1_r,min_X_img2_r);
max_X_img_r = max(max_X_img1_r,max_X_img2_r);


% Interpolation der Bilder an der richtigen Stelle/Pixel - zwei Bilder in 1
% Bild
img1_r = zeros(length(min_Y_img_r:max_Y_img_r),length(min_X_img_r:max_X_img_r));
img2_r = img1_r;
img1_r((min(corners_img1_r(2,:)):max(corners_img1_r(2,:)))+1-min_Y_img_r,...
       (min(corners_img1_r(1,:)):max(corners_img1_r(1,:)))+1-min_X_img_r) ...
       = interp2(X_img1,Y_img1,double(img1),X_img1_dr,Y_img1_dr,'linear',NaN);
img2_r((min(corners_img2_r(2,:)):max(corners_img2_r(2,:)))+1-min_Y_img_r,...
       (min(corners_img2_r(1,:)):max(corners_img2_r(1,:)))+1-min_X_img_r) ...
       = interp2(X_img2,Y_img2,double(img2),X_img2_dr,Y_img2_dr,'linear',NaN);
   
img_r = max(img1_r,img2_r);

% Zwei einzelne Bilder
img_1u_r = zeros(length(min_Y_img1_r:max_Y_img1_r),length(min_X_img1_r:max_X_img1_r));
img_1u_r((min_Y_img1_r:max_Y_img1_r)+1-min_Y_img1_r,...
       (min_X_img1_r:max_X_img1_r)+1-min_X_img1_r) ...
       = interp2(X_img1,Y_img1,double(img1),X_img1_dr,Y_img1_dr,'linear',NaN);
img_2u_r = zeros(length(min_Y_img2_r:max_Y_img2_r),length(min_X_img2_r:max_X_img2_r));
img_2u_r((min_Y_img2_r:max_Y_img2_r)+1-min_Y_img2_r,...
       (min_X_img2_r:max_X_img2_r)+1-min_X_img2_r) ...
       = interp2(X_img2,Y_img2,double(img2),X_img2_dr,Y_img2_dr,'linear',NaN);
figure;
% subplot(1,2,1)
im1 = imshow(uint8(img_r)); hold all;
% subplot(1,2,2)
% im2 = imshow(uint8(img2_r));
% im1.AlphaData = 0.5;
% im2.AlphaData = 0.5;


% Einzelne Bilder übereinander legen
[size_y_im1, size_x_im1] = size(img_1u_r);
[size_y_im2, size_x_im2] = size(img_2u_r);
diff_y = abs(size_y_im1-size_y_im2);
diff_x = abs(size_x_im1-size_x_im2);


figure
subplot(1,2,1)
imshow(uint8(img_1u_r));
subplot(1,2,2)
imshow(uint8(img_2u_r));

figure
im1 = imshow(uint8(img_1u_r));hold all;
im2 = imshow(uint8(img_2u_r));
im1.AlphaData = 0.5;
im2.AlphaData = 0.5;

end
