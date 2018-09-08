function [new_rectified_img,ScatInterp_init] = new_rectified_view(I1_rec,depth_maps,T_x,min_x,max_x,min_y,max_y,cut,p,ScatInterp)
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


if isempty(ScatInterp)
    for i = 1:2
        [OLD_pixel_x,OLD_pixel_y] = meshgrid(1:size(I1_rec{i},2),1:size(I1_rec{i},1));
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
        NEW_pixel_x{i} = OLD_pixel_x + T_x(i)./depth_maps{i};
        % T_y = 0 und T_z = 0, daraus folgt, dass die y- und z-Pixelkoordinaten
        % unverändert bleiben. (z-Komponente hat sowohl in homogenen
        % Bildkoordinaten als auch homogenen Pixelkoordinaten den Wert 1.)
        
        if i == 2
            % Justieren der Pixelkoordinaten => Pixelkoordinaten von rechtem
            % Ausgangsbild auf jene vom linken Ausgangsbild umrechnen.
            NEW_pixel_x{2} = NEW_pixel_x{2}+min_x(1)-min_x(2);
        end
        
        % Valide Pixel bestimmen
        p_valid{i} = isfinite(NEW_pixel_x{i}) & isfinite(depth_maps{i})...
            & NEW_pixel_x{i} >= min_x(1) & NEW_pixel_x{i} <= max_x(1);
        NEW_pixel_x{i} = round(NEW_pixel_x{i}(p_valid{i}));
        NEW_pixel_y{i} = round(OLD_pixel_y(p_valid{i}));
        v{i} = double(I1_rec{i}(p_valid{i}));
        
        if i == 2
            % Justieren der Pixelkoordinaten => Pixelkoordinaten von rechtem
            % Ausgangsbild auf jene vom linken Ausgangsbild umrechnen.
            NEW_pixel_y{2} = round(NEW_pixel_y{2}-min_y(1)+min_y(2));
        end
        
        
        % Doppelte Einträge (Einträge mit gleichen Pixelkoordinaten wegen
        % ungenauer Depth Map) jeweils mitteln.
        hol_grid{i} = [NEW_pixel_x{i},NEW_pixel_y{i}];
        % Bestimme Reihenindizes der validen Pixelkoordinaten (ohne Dopplungen)
        [U_pixel_xy,idx_firstU{i},~] = unique(hol_grid{i},'rows');
        % Mittelwerte bilden
        % accumarray Pixelkoordinaten übergeben mit Ursprung (1,1)
        konst_x{i} = -min(NEW_pixel_x{i})+1;
        konst_y{i} = -min(NEW_pixel_y{i})+1;
        mean_grid = accumarray([hol_grid{i}(:,1)+konst_x{i},hol_grid{i}(:,2)+konst_y{i}],v{i},[],@mean);
        NEW_pixel_x{i} = U_pixel_xy(:,1);
        NEW_pixel_y{i} = U_pixel_xy(:,2);
        v{i} = mean_grid(sub2ind(size(mean_grid),...
            hol_grid{i}(idx_firstU{i},1)+konst_x{i},hol_grid{i}(idx_firstU{i},2)+konst_y{i}));
    end
    
    
    % Sind Informationen aus beiden Bildern für das selbe Pixel erhältlich werden die
    % entsprechenden Intensitäten mit p bzw. (1-p) gewichtet, sodass die
    % Intensitäten der jeweils nähere Ansicht stärker gewichtet werden.
    [log_twice,idx_2] = ismember([NEW_pixel_x{1},NEW_pixel_y{1}],...
        [NEW_pixel_x{2},NEW_pixel_y{2}],'rows');
    for k = 1:length(idx_2)
        % Gewichtung mit p, wenn Pixelkoordinate in beiden Bildern
        % vorhanden.
        if idx_2(k) ~= 0 
            v{2}(idx_2(k)) = p*v{2}(idx_2(k))+(1-p)*v{1}(k);
        end
    end
           
    
    % Entsprechende Daten in v{1} löschen, da nun in v{2}
    % enthalten/berücksichtigt.
    NEW_pixel_x{1}(log_twice) = [];
    NEW_pixel_y{1}(log_twice) = [];
    v{1}(log_twice) = [];
    
    
    % Verwende Informationen aus beiden Bilder, sodass Bereiche, die aus
    % der linken Ansicht nicht sichtbar sind durch die rechte Ansicht
    % bestimmt werden und anders herum.
    NEW_pixel_x = [NEW_pixel_x{1};NEW_pixel_x{2}];
    NEW_pixel_y = [NEW_pixel_y{1};NEW_pixel_y{2}];
    v = [v{1};v{2}];
    
    disp('Creating Scattered Interpolant');

     F = scatteredInterpolant(NEW_pixel_x,NEW_pixel_y,v,'linear','none');
     new_rectified_img = uint8(F({min_x(1):max_x(1),...
         -cut:max_y(1)-min_y(1)-cut}))';

     disp('Scattered Interpolant created');
    
    % Daten speichern um nächsten Farbkanal schneller berechnen zu können.
    for i = 1:2
    ScatInterp_init{i} = struct;
    ScatInterp_init{i}.p_valid = p_valid{i};
    ScatInterp_init{i}.hol_grid = hol_grid{i};
    ScatInterp_init{i}.konst_x = konst_x{i};
    ScatInterp_init{i}.konst_y = konst_y{i};
    ScatInterp_init{i}.idx_firstU = idx_firstU{i};
    end
    ScatInterp_init{1}.idx_2 = idx_2;
    ScatInterp_init{1}.log_twice = log_twice;
    ScatInterp_init{3} = F;
else
    % Scattered Interpolant wurde bereits erstellt (für anderen Farbkanal)
    % => Berechne nur die Intensitäten neu
    for i = 1:2
        % wie oben
        v{i} = double(I1_rec{i}(ScatInterp{i}.p_valid));
        mean_grid = accumarray([ScatInterp{i}.hol_grid(:,1)+ScatInterp{i}.konst_x,...
            ScatInterp{i}.hol_grid(:,2)+ScatInterp{i}.konst_y],v{i},[],@mean);
        v{i} = mean_grid(sub2ind(size(mean_grid),...
           ScatInterp{i}.hol_grid(ScatInterp{i}.idx_firstU,1)+ScatInterp{i}.konst_x,...
           ScatInterp{i}.hol_grid(ScatInterp{i}.idx_firstU,2)+ScatInterp{i}.konst_y));
    end
        for k = 1:length(ScatInterp{1}.idx_2)
            % Gewichtung mit p, wenn Pixelkoordinate in beiden Bildern
            % vorhanden.
            if ScatInterp{1}.idx_2(k) ~= 0 
                v{2}(ScatInterp{1}.idx_2(k)) = p*v{2}(ScatInterp{1}.idx_2(k))+(1-p)*v{1}(k);
            end
        end
        v{1}(ScatInterp{1}.log_twice) = [];
        
        F = ScatInterp{3};
        F.Values = [v{1};v{2}];

        disp('Changing values of Scattered Interpolant');
        
        new_rectified_img = uint8(F({min_x(1):max_x(1),...
             -cut:max_y(1)-min_y(1)-cut}))';

         disp('Values of Scattered Interpolant changed');
        ScatInterp_init = {};
    

end

