function new_rectified_img = new_rectified_view(img_L_rec,img_R_rec,...
    depth_map_L,depth_map_R,T_x,min_L_x,max_L_x,min_R_x,max_R_x)
%NEW_RECTIFIED_VIEW Gibt rektifizierte Zwischenansicht zurück
% -----------------------------------------------------------------------
% !!! T_x zwischen 0 (linke Ansicht) und -1 (rechte Ansicht) wählen !!!
% -----------------------------------------------------------------------
% min_x und max_x bestimmen die Pixelgrenzen für die die neue Ansicht
% berechnet werden soll.

% P2 = R*P1+T mit R=I und P1=lamda1*x1
% lambda2*x2 = lambda1*x1+T => lambda1 und lambda2 sind hier in diesem Fall
% identisch. (Verwendung von rektifizierten Bildern bzw. eukl. Bewegung
% lässt z-Komponente unverändert, x1 und x2 auf z=1 normiert)
% Daher gilt: lambda1 = lambda2 = z-Komponente der Punkte im Raum
% Letzeres ist gegeben durch die Tiefenkarte/depth_map (Achtung Tiefenkarte
% ist in Pixeleinheiten)

% Analoge Rechnung für Pixelkoordinaten:
% => lambda2*x2_pixel = lambda1*K*R*inv(K)*x1_pixel+K*T
% mit K*R*inv(K) = I und K*T = T_pixel
% => lambda2*x2_pixel = lambda1*x1_pixel+T_pixel <= Gleichung(1)

% Nur T_x nötig, da nur Verschiebung in x-Richtung (T = [T_x; 0; 0])
% T_x soll zwischen 0 und -1 liegen

[L_pixel_x,L_pixel_y] = meshgrid(1:size(img_L_rec,2),1:size(img_L_rec,1));

% Tiefe wurde bisher bestimmt als: depth = 1/(x1_pixel-x2_pixel). Der erste
% Eintrag der Kalibrierungsmatrix K(1,1) entspricht dem Umrechnungsfaktor
% zwischen Abständen der x-Weltkoordinaten und x-Pixelkoordinaten:
% d_x_pixel = K(1,1)*d_x_welt
% Daraus folgt für die "echte" Tiefe: depth = K(1,1)/(x1_pixel-x2_pixel)
% Daraus folgt:
% lambda1 = lambda2 = lambda = z-Komponente der Punkte im Raum
%         = K(1,1)*depth_map'
% Umformen von Gleichung(1) ergibt die neuen Pixelkoordinaten:
% x2_pixel = x1_pixel + T_pixel/lambda
% => x2_pixel = x1_pixel + (K(1,1)*T_x) / ((K(1,1))*depth_map)
%             = x1_pixel + T_x / depth_map
NEW_pixel_x = L_pixel_x + T_x./depth_map_L;
% T_y = 0 und T_z = 0, daraus folgt, dass die y- und z-Pixelkoordinaten
% unverändert bleiben. (z-Komponente hat sowohl in homogenen
% Bildkoordinaten als auch homogenen Pixelkoordinaten den Wert 1.)

% Mit der Funktion map_pixel werden in jeder Zeile die Pixel aus den
% rektifizierten Bildern in die virtuelle Ansicht gemapped
% img_map = map_pixel(img_L_rec,[],depth_map_L,[],T_x,[],min_L_x,max_L_x,[],[]);

% Valide Pixel bestimmen
p_valid = isfinite(NEW_pixel_x) & isfinite(depth_map_L);
NEW_pixel_x = NEW_pixel_x(p_valid);
NEW_pixel_y = L_pixel_y(p_valid);
v = double(img_L_rec(p_valid));

F = scatteredInterpolant(NEW_pixel_x,NEW_pixel_y,v,'linear','none');
new_rectified_img = uint8(F({min_L_x:max_L_x,1:size(img_L_rec,1)}))';
end

