function [virtual_view_img,img_rectified_new] = projection(IGray_L,...
                        img_rectified_L_full,depth_map,K,T_x,R_rect,...
                        alpha,offset_x,d_cut_up,d_cut_down)
%PROJECTION Berechnet neue rektifizierte Zwischenansicht und führt
%anschließend Derektifizierung durch
% IGray_L = Linkes Graubild
% img_rectified_L = Linkes rektifiziertes Bild
% depth_map = Tiefenkarte (in Pixeleinheiten)
% K = Kalibrierungsmatrix
% T_x zwischen 0 (linke Ansicht) und -1 (rechte Ansicht) wählen. (siehe
% new_rectified_view.m)
% R_rect = Zuvor berechnete Rotationsmatrix für die Rektifizierung
% min_x_pixel_rec = Minimale unkalibrierte x-Weltkoordinate der
% rektifizierten Zwischenansicht. (siehe new_rectified_view.m)
% max_x_pixel_rec = Maximale unkalibrierte x-Weltkoordinate der
% rektifizierten Zwischenansicht. (siehe new_rectified_view.m)
% alpha = Winkel, um den die Ansicht vom linken Bild ausgehend um die
% y-Achse gedreht wurde.

% Gewünschte Rotation um y-Achse bestimmen
R_v = [cosd(alpha),0,sind(alpha);
       0,1,0;
       -sind(alpha),0,cosd(alpha)];
% Umwandlung der derektifizierten Pixelkoordinaten x in rektifizierte
% Pixelkoordinaten x':
% x' = K*R_rect'*R_v*inv(K)*x
% ("Derektifizierte Ansicht zunächst um R_rect nach rechts drehen und
% anschließend mit R_v um Winkel alpha nach links drehen")
% Bzw. inverse Transformation:
% Umwandlung der rektifizierten Pixelkoordinaten x' in derektifizierte
% Pixelkoordinaten x:
% x = K*R_rect*R_v'*inv(K)*x
% ("Rektifizierte Ansicht zunächst um R_rect nach links drehen und
% anschließend mit R_v um Winkel alpha nach rechts drehen")
T_rec = K*R_rect*R_v'*inv(K);

% Berechnung & Drehung der Bildmitte, so dass Bildebenen parallel zu
% Translationsvektor liegen
[size_y, size_x] = size(IGray_L);
center_y = round(size_y/2);
center_x = round(size_x/2);
center = [center_x; center_y; 1];
center_rec = T_rec*center;

% Horizontaler/Vertikaler Versatz der begradigten Bilder soll nun gleich
% sein, bzw. Bilder sollen horizontal/vertikal auf gleicher Höhe liegen.
d_rec = center(1:2) - center_rec(1:2)/center_rec(3);
% Anpassen der Kalibrierungsmatrix
K_rec = K;
% Anpassen des Ursprungs in x-Richtung
K_rec(1,3) = K(1,3)+d_rec(1);
% Anpassen des Ursprungs in y-Richtung
K_rec(2,3) = K(2,3)+d_rec(2);
% Neue Transformation für die Rektifizierung berechnen.
T_rec = K_rec*R_rect*R_v'*inv(K);

% Bestimme wo die Ecken und Mittelpunkt der derektifizierten Ansicht in der
% rektifizierten Ansicht liegen. 
corners_rec = T_rec*[0, size_x,size_x, 0;
                size_y, size_y, 0, 0;
                ones(1,4)];
center_rec = T_rec*center;
% Normalisierung der z-Komponente auf eins.
center_rec = center_rec/center_rec(3);
corners_rec = corners_rec./corners_rec(3,:);

% Bestimme minimale und maximale x-Koordinate von corners_rec
min_x_v_rec = min(corners_rec(1,:));
max_x_v_rec = max(corners_rec(1,:));

% Ansicht wird nach rechts verschoben. Diese Verschiebung nach rechts
% bedeutet, dass der neue Pixelbereich, für den die neuen Intensitäten zu
% berechnen sind, nach links verschoben wird.
% Variable dif entspricht dem Abstand in x-Richtung zwischen der linken
% Ecke der neuen rektifizierten Ansicht und ihrem ursprünglichen
% Mittelpunkt (vor der Vorschiebung um d_rec(1)).
size_x_rec = size(img_rectified_L_full,2);
dif = center_rec(1)-min_x_v_rec+d_rec(1);
% Sicher gehen, dass min_x_v_rec <= 0
min_x_v_rec = min(-round(dif-(size_x_rec-offset_x)/2),0);
max_x_v_rec = max_x_v_rec+min_x_v_rec;


% Bestimme neue rektifizierte Zwischenansicht
img_rectified_new = new_rectified_view(img_rectified_L_full,depth_map,T_x,...
    min_x_v_rec,max_x_v_rec);

% Meshgrid, das Pixelkoordinaten der derektifizierten Ansicht enthält in
% rektifizierte Pixelkoodrinaten transformieren.
[pixel_x, pixel_y] = meshgrid(1:size_x,1:size_y);
% Entsprechende rektifizierte Pixelkoordinaten berechnen.
% Normierung auf z = 1 => Dazu als erstes z-Komponente berechnen.
z_rec = T_rec(3,1)*pixel_x+T_rec(3,2)*pixel_y+T_rec(3,3);
pixel_x_rec = (T_rec(1,1)*pixel_x+T_rec(1,2)*pixel_y+T_rec(1,3)) ./ z_rec;
pixel_y_rec = (T_rec(2,1)*pixel_x+T_rec(2,2)*pixel_y+T_rec(2,3)) ./ z_rec;

% Interpolation
v = double(img_rectified_new);
% Um Interpolant zu bestimmen, bekannte Stützstellen und Intensitäten
% übergeben.
F = griddedInterpolant({1:size(img_rectified_new,2),...
        1-d_cut_up:size(img_rectified_new,1)-d_cut_up},v','linear','none');
virtual_view_img = uint8(F(pixel_x_rec(:),pixel_y_rec(:)))';
virtual_view_img = reshape(virtual_view_img,[size_y, size_x]);
virtual_view_img = virtual_view_img(1:end-d_cut_down,:);
end
