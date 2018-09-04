function [virtual_image,virtual_image_rectified] = projection_b(img_L,img_R,...
    img_L_rect_full,img_R_rect_full,...
    depth_map_L,depth_map_R,K,T,R,p)
%PROJECTION_B Summary of this function goes here
%   Detailed explanation goes here

%% Bestimme Rotation und Translation f�r Zwischenansicht
[theta_deg, psi_deg, phi_deg] = euler_angles(R);
theta_deg_fv =  p*theta_deg;
psi_deg_fv  = p*psi_deg;
phi_deg_fv = p*phi_deg;
R_v = euler_rotation(phi_deg_fv, theta_deg_fv, psi_deg_fv);
T_v = p*T;

%% Rectification Algorithum (von Paper Fusiello)
% Mit der zwischen Ansicht sind 3 Kamerazentren vorhanden PoFV ist nun die
% Zwischenansicht zwischen Po1 und Po2
Po1 = K*[1 0 0 0; 0 1 0 0; 0 0 1 0];
PoFV = K*[R_v, T_v];
Po2 = K*[R, T];
[TL_L_FV, TR_L_FV, R_rect_L_FV] = rectify_fusiello(Po1,PoFV,K);
[TL_FV_R, TR_FV_R, R_rect_FV_R] = rectify_fusiello(PoFV,Po2,K);

% Schnittpunkt optische Achse mit Bildebene entspricht Bildmitte
% Linkes und rechtes Bild sind gleich gross.
[row, col] = size(img_L);
center_row = round(row/2);
center_col = round(col/2);

% Berechnung & Drehung der Bildmitte, so dass Bildebenen parallel zu
% Translationsvektor liegen
center = [center_col; center_row; 1];
center_L_r = TL_L_FV*center;
center_L_FV_r = TR_L_FV*center;
center_R_FV_r = TL_FV_R*center;
center_R_r = TR_FV_R*center;


% Vertikaler Versatz der begradigten Bilder (Y-Achse) soll nun gleich
% sein, bzw. Bilder sollen vertikal auf gleicher Hoehe liegen.
dL(2) = center(2) - center_L_r(2)/center_L_r(3);
dL_FV(2) = center(2) - center_L_FV_r(2)/center_L_FV_r(3);
dR_FV(2) = center(2) - center_R_FV_r(2)/center_R_FV_r(3);
dR(2) = center(2) - center_R_r(2)/center_R_r(3);


% Bildausschnitt in x-Richtung wird spaeter passend gewaehlt.
dL(1) = 0;
dL_FV(1) = 0;
dR_FV(1) = 0;
dR(1) = 0;

% Homographien erneut mit angepasster Kalibrierung berechnen.
[TL_L_FV, TR_L_FV, R_rect_L_FV] = rectify_fusiello(Po1,PoFV,K,R_rect_L_FV,dL,dL_FV);
[TL_FV_R, TR_FV_R, R_rect_FV_R] = rectify_fusiello(PoFV,Po2,K,R_rect_FV_R,dR_FV,dR);

TL_L_FV_inv = inv(TL_L_FV);
TR_L_FV_inv = inv(TR_L_FV);
TL_FV_R_inv = inv(TL_FV_R);
TR_FV_R_inv = inv(TR_FV_R);


% Zur Interpolation werden die Bildaussschnitte zuerst in das neue
% Koordinatensystem und daraufhin zurücktransformiert. Hierfür müssen
% die Bildecken bekannt sein.
% Linkes und rechtes Bild sind gleich groß.
[size_y, size_x] = size(img_L);
center_L_r = TL_L_FV*center;
center_L_FV_r = TR_L_FV*center;
center_R_FV_r = TL_FV_R*center;
center_R_r = TR_FV_R*center;
% Bildecken in rotierten Pixelkoordinaten
corners_L_r = TL_L_FV*[0, size_x,size_x, 0;
    size_y, size_y, 0, 0;
    ones(1,4)]; 
corners_L_r = corners_L_r./corners_L_r(3,:);% Normalisierung der z-Komponente auf eins.

center_L_FV_r = TR_L_FV*[0, size_x,size_x, 0;
    size_y, size_y, 0, 0;
    ones(1,4)]; 
center_L_FV_r = center_L_FV_r./center_L_FV_r(3,:);% Normalisierung der z-Komponente auf eins.

center_R_FV_r = TL_FV_R*[0, size_x,size_x, 0;
    size_y, size_y, 0, 0;
    ones(1,4)]; 
center_R_FV_r = center_R_FV_r./center_R_FV_r(3,:);% Normalisierung der z-Komponente auf eins.

corners_R_r = TR_FV_R*[0, size_x,size_x,0;
    size_y, size_y,0,0;
    ones(1,4)];
corners_R_r = corners_R_r./corners_R_r(3,:);% Normalisierung der z-Komponente auf eins.

min_x_L_FV = ceil(min(center_L_FV_r(1,:)));
max_x_L_FV = floor(max(center_L_FV_r(1,:)));
min_x_R_FV = ceil(min(center_R_FV_r(1,:)));
max_x_R_FV = floor(max(center_R_FV_r(1,:)));

[img_map_L,img_map_R] = map_pixel(img_L_rect_full,img_R_rect_full,...
                    depth_map_L,depth_map_R,p,min_x_L_FV,max_x_L_FV,...
                    min_x_R_FV,max_x_R_FV);


end

