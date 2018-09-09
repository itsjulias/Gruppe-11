function [disparity_map,hom_map] = disparity_estimation(I_L_rec,I_R_rec,...
    interv_search_left,interv_search_right,off_d,window_length,L_R_base)
%DEPTH_ESTIMATION2 mittels SAD Liniensuche

% I_L_rec und I_R_rec haben selbe Anzahl an Zeilen, aber nicht unbedingt selbe
% Anzahl an Spalten
[zl_L,sp_L] = size(I_L_rec);
sp_R = size(I_R_rec,2);

disparity_map = zeros(sp_L,zl_L);
% hom_map: Enthält Nullen an Stellen, an denen eine homogene Fläche
% detektiert wurde, sonst Einsen.
hom_map = zeros(sp_L,zl_L);

% Erlaubte minimale Distanz zum Bildrand.
s = (window_length-1)/2;

% Liniensuche zeilenweise
for i = s+1:zl_L-s % Zeilen werden hochgezÃ¤hlt, s: Abstand zum Rand, zl_L: unterer Rand Bild1
    [disparity_map(s+1:sp_L-s,i),hom_map(s+1:sp_L-s,i)] = SAD_line(i,L_R_base); %sp_L: rechter Rand Bild1
end

% Ränder füllen und transponieren
% Linken Rand füllen
disparity_map(:,1:s) = repmat(disparity_map(:,s+1),[1 s]);
% Rechten Rand füllen
disparity_map(:,end-s+1:end) = repmat(disparity_map(:,end-s),[1 s]);
% Oberen Rand füllen
disparity_map(1:s,:) = repmat(disparity_map(s+1,:),[s 1]);
% Unteren Rand füllen
disparity_map(end-s+1:end,:) = repmat(disparity_map(end-s,:),[s 1]);
% Rückgabe der Disparity Map entsprechend des linken rektifizierten Bildformats
disparity_map = disparity_map';
hom_map = hom_map';

function [d_cor,line_hom] = SAD_line(line_number,L_R_base) % , min_cor
        % Function gibt nur d fÃƒÂ¯Ã‚Â¿Ã‚Â½r die Pixel des linkes Bildes mit mind.
        % Abstand s zum Bildrand zurÃƒÂ¯Ã‚Â¿Ã‚Â½ck => length(d_cor) = sp-2*s
        % log_hom: Enthält Nullen an Stellen, an denen eine homogene Fläche
        % detektiert wurde, sonst Einsen.
        line_hom = ones(1,sp_L-2*s);
        
        % Matrizen fÃ¼r stacked vectoren anlegen.
        N_L = zeros(window_length^2,sp_L-2*s);
        N_R = zeros(window_length^2,sp_R-2*s);
        % Normierte Fenster fÃƒÂ¼r linkes Bild berechnen:
        for k = s+1:sp_L-s
             % Quadratisches Fenster extrahieren.
             L_W = double(I_L_rec(line_number-s:line_number+s,k-s:k+s));
             idx = k-s;
             if std(L_W(:)) < 3
                line_hom(idx) = 0; 
             end
             N_L(:,idx) = L_W(:);
        end
        % analog fÃƒÂ¼r rechtes Bild
        for k = s+1:sp_R-s
             R_W = double(I_R_rec(line_number-s:...
                 line_number+s,k-s:k+s));
             idx = k-s;
             N_R(:,idx) = R_W(:);
        end
        % Suche fÃƒÂ¯Ã‚Â¿Ã‚Â½r jedes Pixel im linken Bild passendes Pixel im rechten
        % Bild in einem Intervall von max_d bis min_d Pixel.
        % NCC-matrix: Zeile entspricht x-Pixelkoordinate im linken Bild,
        % Spalte entspricht dem Liniensuchintervall im rechten Bild.
        SAD_matrix = zeros(sp_L-2*s,interv_search_right-interv_search_left+1);

        % Iterieren durch alle Spalten von N_L bzw. N_R
        for X = 1:size(N_L,2)
            N_L_W = N_L(:,X);
            for Y = interv_search_left:interv_search_right
                if X+Y < 1 || X+Y > size(N_R,2)
                    % Pixel in Bild 2 liegt nicht innerhalb des Bildes mit
                    % mindestens Abstand s zum Bildrand
                    SAD_matrix(X,Y-interv_search_left+1) = Inf;
                else
                    % Korrelation mittels Skalarprodukt berechnen.
                    SAD_matrix(X,Y-interv_search_left+1) = sum(abs(N_L_W-N_R(:,X+Y)));
                end
            end
        end
        [~, ind_min] = min(SAD_matrix,[],2);   

        % Spaltenvektor erstellen, der die PixelabstÃƒÂ¯Ã‚Â¿Ã‚Â½nde zwischen einem
        % Pixel im linken Bild und dem korrelierten Pixel im rechten Bild
        % enthÃƒÂ¯Ã‚Â¿Ã‚Â½lt. Erster Eintrag enthÃƒÂ¯Ã‚Â¿Ã‚Â½lt d vom Pixel mit der x-Koordinate
        % s+1 in der aktuellen Suchzeile im linken Bild.
        % d = m1 - m2 (Offset berÃƒÂ¼cksichtigen!)
        if strcmp(L_R_base,'R_base')
            d_cor = off_d+interv_search_left+ind_min-1;
        else
            d_cor = off_d-(ind_min+interv_search_left-1);
        end
    end



function d_cor = NCC_line(line_number) % , min_cor
        % Function gibt nur d fÃƒÂ¯Ã‚Â¿Ã‚Â½r die Pixel des linkes Bildes mit mind.
        % Abstand s zum Bildrand zurÃƒÂ¯Ã‚Â¿Ã‚Â½ck => length(d_cor) = sp-2*s
               
        % Matrizen fÃƒÂ¼r stacked vectoren anlegen.
        N_L = zeros(window_length^2,sp_L-2*s);
        N_R = zeros(window_length^2,sp_R-2*s);
        % Normierte Fenster fÃƒÂ¼r linkes Bild berechnen:
        for k = s+1:sp_L-s
             % Quadratisches Fenster extrahieren.
             L_W = double(I_L_rec(line_number-s:line_number+s,k-s:k+s));
             % Mittelwert subtrahieren.
             L_W = L_W - mean(L_W(:));
             % Teilen durch die Standardabweichung
             L_W = L_W/std(L_W(:));
             % Normierte Matrix als Vektor speichern.
             idx = k-s;
             N_L(:,idx) = L_W(:);
        end
        % analog fÃƒÂ¼r rechtes Bild
        for k = s+1:sp_R-s
            % Suche im rechten Bild eine Zeile ÃƒÂ¼ber der Zeile des linken
            % Bilds
             R_W = double(I_R_rec(line_number-s:...
                 line_number+s,k-s:k+s));
             R_W = R_W - mean(R_W(:));
             R_W = R_W/std(R_W(:));
             idx = k-s;
             N_R(:,idx) = R_W(:);
        end
        % Suche fÃƒÂ¯Ã‚Â¿Ã‚Â½r jedes Pixel im linken Bild passendes Pixel im rechten
        % Bild in einem Intervall von max_d bis min_d Pixel.
        % NCC-matrix: Zeile entspricht x-Pixelkoordinate im linken Bild,
        % Spalte entspricht dem Liniensuchintervall im rechten Bild.
        NCC_matrix = zeros(sp_L-2*s,interv_search_right-interv_search_left+1);

        % Iterieren durch alle Spalten von N_L bzw. N_R
        for X = 1:size(N_L,2)
            N_L_W = N_L(:,X);
            for Y = interv_search_left:interv_search_right
                if X+Y < 1 || X+Y > size(N_R,2)
                    % Pixel in Bild 2 liegt nicht innerhalb des Bildes mit
                    % mindestens Abstand s zum Bildrand
                    NCC_matrix(X,Y-interv_search_left+1) = -1;
                else
                    % Korrelation mittels Skalarprodukt berechnen.
                    NCC_matrix(X,Y-interv_search_left+1) = (1/(window_length^2-1))*(N_L_W'*N_R(:,X+Y));
                end
            end
        end
        [~, ind_max] = max(NCC_matrix,[],2);   

        % Spaltenvektor erstellen, der die PixelabstÃƒÂ¯Ã‚Â¿Ã‚Â½nde zwischen einem
        % Pixel im linken Bild und dem korrelierten Pixel im rechten Bild
        % enthÃƒÂ¯Ã‚Â¿Ã‚Â½lt. Erster Eintrag enthÃƒÂ¯Ã‚Â¿Ã‚Â½lt d vom Pixel mit der x-Koordinate
        % s+1 in der aktuellen Suchzeile im linken Bild.
        % d = m1 - m2 (Offset berÃƒÂ¼cksichtigen!)
        d_cor = off_d-(ind_max+interv_search_left-1);
    end

end
