function [virtual_image] = projection_b(img_L,img_R,...
    img_L_rect_full,img_R_rect_full,...
    depth_map_L,depth_map_R,K,T,R,p,min_y_full,max_y_full)
%PROJECTION_B Projiziert die beiden rektifizierten Ansichten
%(links,rechts) in die virtuelle Ansicht. Die depth_maps müssen dabei
%sowohl für das linke, als auch das rechte Bild positive Werte beinhalten.
%Die beiden Bilder img_L (und img_R) dienen nur zur Bestimmung der Größe der
%virutellen Ansicht. min_y_fiull, max_y_full sind jeweils die kleinste und
%die größte y-Koordinate der Bilder img_L_rect_full und img_R_rect full
%(beide rektifizierten Bilder müssen in y-Richtung gleich groß sein).

%% Bestimme Rotation und Translation für Zwischenansicht
% Um eine anteilige Rotation und Translation für die Zwischenansicht
% generieren zu können, müssen zunächst die Euler-Winkel aus der Rotation R
% bestimmt werden um dieser darauf folgend zu skalieren und in die
% Rotationsmatrix der virtuellen Ansicht umzurechnen.
[theta_deg, psi_deg, phi_deg] = euler_angles(R);
theta_deg_fv =  p*theta_deg;
psi_deg_fv  = p*psi_deg;
phi_deg_fv = p*phi_deg;
R_v = euler_rotation(phi_deg_fv, theta_deg_fv, psi_deg_fv); % Rotation der virtuellen Ansicht
T_v = p*T; %Translation der virtuellen Ansicht

%% Rectification Algorithum (von Paper Fusiello)
% Mit der virtuellen Ansicht, der linken Ansicht und der rechten Ansicht
% sind 3 Kamerazentren vorhanden. PoFV ist nun die Zwischenansicht zwischen
% Po1 (links) und Po2 (rechts)
Po1 = K*[1 0 0 0; 0 1 0 0; 0 0 1 0]; %linke Kamera: Rotation = Einheitsmatrix, Translation = 0
PoFV = K*[R_v, T_v]; %virtuelle Kamera: Rotation = R_v, Translation =  T_v
Po2 = K*[R, T];%rechte Kamera: Rotation = R, Translation = T

%Berechnung der Rektifizierungs-Transformationen
[TL_L_FV, TR_L_FV, R_rect_L_FV] = rectify_fusiello(Po1,PoFV,K);
[TL_FV_R, TR_FV_R, R_rect_FV_R] = rectify_fusiello(PoFV,Po2,K);

% Schnittpunkt optische Achse mit Bildebene entspricht Bildmitte
% Linkes und rechtes Bild sind gleich gross. (s. rectification.m)
[row, col] = size(img_L);
center_row = round(row/2);
center_col = round(col/2);

% Berechnung & Drehung der Bildmitte, so dass Bildebenen parallel zu
% Translationsvektor liegen (s. rectification.m)
center = [center_col; center_row; 1];
center_L_r = TL_L_FV*center;
center_L_FV_r = TR_L_FV*center;
center_R_FV_r = TL_FV_R*center;
center_R_r = TR_FV_R*center;


% Vertikaler Versatz der begradigten Bilder (Y-Achse) soll nun gleich
% sein, bzw. Bilder sollen vertikal auf gleicher Hoehe liegen. (s.
% rectification.m) 
dL(2) = center(2) - center_L_r(2)/center_L_r(3);
dL_FV(2) = center(2) - center_L_FV_r(2)/center_L_FV_r(3);
dR_FV(2) = center(2) - center_R_FV_r(2)/center_R_FV_r(3);
dR(2) = center(2) - center_R_r(2)/center_R_r(3);


% Bildausschnitt in x-Richtung wird spaeter passend gewaehlt. (s.
% rectification.m) 
dL(1) = 0;
dL_FV(1) = 0;
dR_FV(1) = 0;
dR(1) = 0;

% Homographien erneut mit angepasster Kalibrierung berechnen. (s.
% rectification.m) 
[TL_L_FV, TR_L_FV, R_rect_L_FV] = rectify_fusiello(Po1,PoFV,K,R_rect_L_FV,dL,dL_FV);
[TL_FV_R, TR_FV_R, R_rect_FV_R] = rectify_fusiello(PoFV,Po2,K,R_rect_FV_R,dR_FV,dR);

TL_L_FV_inv = inv(TL_L_FV);
TR_L_FV_inv = inv(TR_L_FV);
TL_FV_R_inv = inv(TL_FV_R);
TR_FV_R_inv = inv(TR_FV_R);


% Zur Interpolation werden die Bildaussschnitte zuerst in das neue
% Koordinatensystem und daraufhin zuruecktransformiert. Hierfuer muessen
% die Bildecken bekannt sein.
% Linkes und rechtes Bild sind gleich gross.
[size_y, size_x] = size(img_L);
center_L_r = TL_L_FV*center;
center_L_FV_r = TR_L_FV*center;
center_R_FV_r = TL_FV_R*center;
center_R_r = TR_FV_R*center;
% Bildecken in rotierten Pixelkoordinaten

% für linkes Bild
corners_L_r = TL_L_FV*[0, size_x,size_x, 0;
    size_y, size_y, 0, 0;
    ones(1,4)]; 
corners_L_r = corners_L_r./corners_L_r(3,:);% Normalisierung der z-Komponente auf eins.

% für virtuelle Ansicht (aus Sicht links)
corners_L_FV_r = TR_L_FV*[0, size_x,size_x, 0;
    size_y, size_y, 0, 0;
    ones(1,4)]; 
corners_L_FV_r = corners_L_FV_r./corners_L_FV_r(3,:);% Normalisierung der z-Komponente auf eins.

%für virtuelle Ansicht (aus Sicht rechts)
corners_R_FV_r = TL_FV_R*[0, size_x,size_x, 0;
    size_y, size_y, 0, 0;
    ones(1,4)]; 
corners_R_FV_r = corners_R_FV_r./corners_R_FV_r(3,:);% Normalisierung der z-Komponente auf eins.

%für rechtes Bild
corners_R_r = TR_FV_R*[0, size_x,size_x,0;
    size_y, size_y,0,0;
    ones(1,4)];
corners_R_r = corners_R_r./corners_R_r(3,:);% Normalisierung der z-Komponente auf eins.

% Festlegen der Bildgrenzen (Pixel-Koordinaten in x-Richtung) der virtuellen
% Ansicht (bezogen auf das linke Bild), um später mit der Tiefenkarte die
% Pixel an die richtigen Stellen projizieren zu können.
min_x_L_FV = ceil(min(corners_L_FV_r(1,:)));
max_x_L_FV = floor(max(corners_L_FV_r(1,:)));
% Festlegen in y-Richtung (wird nicht benötigt, da min_y_full und
% max_y_full vorhanden sind)
% min_y_L_FV = ceil(min(corners_L_FV_r(2,:)));
% max_y_L_FV = floor(max(corners_L_FV_r(2,:)));

% Festlegen der Bildgrenzen (Pixel-Koordinaten in x-Richtung) der virtuellen
% Ansicht (bezogen auf das rechte Bild), um später mit der Tiefenkarte die
% Pixel an die richtigen Stellen projizieren zu können.
min_x_R_FV = ceil(min(corners_R_FV_r(1,:)));
max_x_R_FV = floor(max(corners_R_FV_r(1,:)));
% Festlegen in y-Richtung (wird nicht benötigt, da min_y_full und
% max_y_full vorhanden sind)
% min_y_R_FV = ceil(min(corners_R_FV_r(2,:)));
% max_y_R_FV = floor(max(corners_R_FV_r(2,:)));

% In map_pixel werden die einzelnen Pixel aus der linke und rechte
% rektifizierten Ansicht in die virtuelle Ansicht projiziert. (nähere
% Beschreibung s. map_pixel.m)
[img_map_L,img_map_R, img_mapped] = map_pixel(img_L_rect_full,img_R_rect_full,...
                    depth_map_L,depth_map_R,p,min_x_L_FV,max_x_L_FV,...
                    min_x_R_FV,max_x_R_FV);
                
%  img_R_r_full = uint8(interp2(meshX,meshY,double(img2),...
%     meshX_R_dr_full,meshY_R_dr_full,'linear',NaN));

% Hier wird Pixelbereich für virtuelle Ansicht festgelegt. Die Pixel
% koordinaten entsprechen dann 1:2000 in y-Richtung und 1:3000 in
% x-Richtung im Fall von Bildern L1.jpg R1.jpg oder L2.jpg R2.jpg.
meshZ_FV_dr = ones(size(img_L));
[meshX_FV_dr,meshY_FV_dr] = meshgrid(1:size(img_L,2),1:size(img_L,1));

% Die Pixelkoordinaten werden nun in die rektifizierte Ansicht projiziert
% um später img_map in diesem Bereich interpolieren zu können. Die
% Rücktransformation findet quasi via interpolation über die rektifizierten
% Koordinaten statt.
meshZ_FV_r = TR_L_FV(3,1)*meshX_FV_dr + TR_L_FV(3,2)*meshY_FV_dr +...
                  TR_L_FV(3,3)*meshZ_FV_dr;              
meshX_FV_r = (TR_L_FV(1,1)*meshX_FV_dr + TR_L_FV(1,2)*meshY_FV_dr +...
                  TR_L_FV(1,3)*meshZ_FV_dr)./meshZ_FV_r; 
meshY_FV_r = (TR_L_FV(2,1)*meshX_FV_dr + TR_L_FV(2,2)*meshY_FV_dr +...
              TR_L_FV(2,3)*meshZ_FV_dr)./meshZ_FV_r; 
          
% Nun werden die Stuetzstellen von img_mapped festgelegt. Diese sind zum
% einen in y-Richtung durch die min_y_full und max_y_full Koordinaten der
% linken und rechten rektifizierten Ansicht bekannt und zum anderen durch
% die rektifizierten Ecken der virtuellen Ansicht
% (corners_L_FV_r,corners_R_FV_r). Hier wird die Transformation immer
% ausgehend vom linken Bild durchgeführt.
[meshX,meshY] = meshgrid(min_x_L_FV:max_x_L_FV,min_y_full:max_y_full); 

% Fuer Interpolation muss mapped Ansicht in double konvertiert werden
% img_map_dbl = double(img_mapped(1:size(img_map_L,1),1:size(img_map_L,2)));
img_map_dbl = double(img_mapped);

useFilling = false;
if(~useFilling)
% ScatteredInterpolant wird verwendet um die lueckenhafen stellen (also alle
% schwarzen Pixel der gemappten Ansicht img_mapped) zu interpolieren,
% soweit moeglich
F = scatteredInterpolant(meshY(img_map_dbl~=0),meshX(img_map_dbl~=0),img_map_dbl(img_map_dbl~=0),'linear','none');
img_free_viewpoint_interp = uint8(F({meshY_FV_r(1:end,1),meshX_FV_r(1,1:end)}));
else
% Alternativ dazu kann die mapped Ansicht (img_mapped) zuerst mit einem
% fuell-Filter an den loechrigen Stellen gefuellt werden, z.B. median und
% darauf hin mit interp2 interpoliert werden (viel schneller als
% scatteredInterpolant, aber Fuell-Filter dauert entsprechend lange)
img_map_dbl = double(fill_disparity_map(double(img_mapped),img_mapped==0,55,1));
img_free_viewpoint_interp = uint8(interp2(meshX,meshY,...
    img_map_dbl,meshX_FV_r,meshY_FV_r));
end

% figure
% imshow(uint8(img_free_viewpoint_interp));

virtual_image = img_free_viewpoint_interp;
end

