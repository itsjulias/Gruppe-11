function [virtual_view_img,img_rectified_new,ScatInterp_init] = projection(IGray_L,...
                        img_rectified_L_full,img_rectified_R_full,depth_maps,K,R,p,R_rect,...
                        offset_x,min_corners_L_y,ScatInterp)

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

[theta_deg, psi_deg, phi_deg] = euler_angles(R);

for i=1:2
    k = p-i+1;
    if i == 1
         T_x(i) = -p;

         % Bestimme Rotation und Translation für Zwischenansicht
         theta_deg_fv =  k*theta_deg;
         psi_deg_fv  = k*psi_deg;
         phi_deg_fv = k*phi_deg;
         R_v = euler_rotation(phi_deg_fv, theta_deg_fv, psi_deg_fv);

        % Umwandlung der derektifizierten Pixelkoordinaten x in rektifizierte
        % Pixelkoordinaten x':
        % x' = K*R_rect'*R_v*inv(K)*x
        % ("Derektifizierte Ansicht zunächst um R_rect nach rechts drehen und
        % anschließend mit R_v um Winkel alpha nach links drehen")
        % Bzw. inverse Transformation:
        % Umwandlung der rektifizierten Pixelkoordinaten x' in derektifizierte
        % Pixelkoordinaten x:
        % x = K*R_rect*R_v'*inv(K)*x'
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
    
        img_rectified_full{i} = img_rectified_L_full;
    else
        T_x(i) = -k;
        img_rectified_full{i} = img_rectified_R_full;
    end

    % Ansicht wird nach rechts verschoben. Diese Verschiebung nach rechts
    % bedeutet, dass der neue Pixelbereich, für den die neuen Intensitäten zu
    % berechnen sind, nach links verschoben wird.
    if i == 1
        min_x_v_rec(i) = round(-p*offset_x);
    else
        min_x_v_rec(i) = round(-k*offset_x);
    end
    max_x_v_rec(i) = min_x_v_rec(i)+size(img_rectified_full{i},2);
    
    off_y = center_rec(2)-size(img_rectified_full{i},1)/2; 
    min_y_v_rec(i) = round(min(corners_rec(2,:))+off_y);
    max_y_v_rec(i) = round(max(corners_rec(2,:))+off_y);
    
    if i == 1
        cut = min_corners_L_y-min(corners_rec(2,:));
    end
end

if nargin < 11
    % Es sind keine ScatteredInterpolants übergeben worden. Sie müssen neu
    % erstellt werden.
    ScatInterp = [];
end
    % Bestimme neue rektifizierte Zwischenansicht
    [img_rectified_new,ScatInterp_init] = new_rectified_view(img_rectified_full,depth_maps,T_x,...
        min_x_v_rec,max_x_v_rec,min_y_v_rec,max_y_v_rec,cut,p,ScatInterp);

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
    min_corners_x = min(corners_rec(1,:));
    min_corners_y = min(corners_rec(2,:));
    F = griddedInterpolant({1+min_corners_x:size(img_rectified_new,2)+min_corners_x,...
            1+min_corners_y:size(img_rectified_new,1)+min_corners_y},v','linear'); % ,'none'
    virtual_view_img = uint8(F(pixel_x_rec(:),pixel_y_rec(:)))';
    virtual_view_img = reshape(virtual_view_img,[size_y, size_x]);
    

end

