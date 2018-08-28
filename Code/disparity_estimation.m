function disparity_map = disparity_estimation(I_L_rec,I_R_rec,...
    interv_search_left,interv_search_right,off_d,window_length)
%DEPTH_ESTIMATION2 mittels NCC Liniensuche

% I_rec und I2_rec haben selbe Anzahl an Zeilen, aber nicht unbedingt selbe
% Anzahl an Spalten
[zl_L,sp_L] = size(I_L_rec);
sp_R = size(I_R_rec,2);

disparity_map = zeros(sp_L,zl_L);

% Es wird nicht jeweils die komplette Zeile nach Korrespondenzen abgesucht,
% sondern nur innerhalb des festgelegten Intervalls.
r_d = -interv_search_right;
l_d = interv_search_left;

% Erlaubte minimale Distanz zum Bildrand.
s = (window_length-1)/2;

% Liniensuche zeilenweise
for i = s+1:zl_L-s % Zeilen werden hochgezählt, s: Abstand zum Rand, zl_L: unterer Rand Bild1
    disparity_map(s+1:sp_L-s,i) = SAD_line(i); %sp_L: rechter Rand Bild1
end

disparity_map = disparity_map';


    function d_cor = NCC_line(line_number) % , min_cor
        % Function gibt nur d fÃ¯Â¿Â½r die Pixel des linkes Bildes mit mind.
        % Abstand s zum Bildrand zurÃ¯Â¿Â½ck => length(d_cor) = sp-2*s
               
        % Matrizen fÃ¼r stacked vectoren anlegen.
        N_L = zeros(window_length^2,sp_L-2*s);
        N_R = zeros(window_length^2,sp_R-2*s);
        % Normierte Fenster fÃ¼r linkes Bild berechnen:
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
        % analog fÃ¼r rechtes Bild
        for k = s+1:sp_R-s
            % Suche im rechten Bild eine Zeile Ã¼ber der Zeile des linken
            % Bilds
             R_W = double(I_R_rec(line_number-s:...
                 line_number+s,k-s:k+s));
             R_W = R_W - mean(R_W(:));
             R_W = R_W/std(R_W(:));
             idx = k-s;
             N_R(:,idx) = R_W(:);
        end
        % Suche fÃ¯Â¿Â½r jedes Pixel im linken Bild passendes Pixel im rechten
        % Bild in einem Intervall von max_d bis min_d Pixel.
        % NCC-matrix: Zeile entspricht x-Pixelkoordinate im linken Bild,
        % Spalte entspricht dem Liniensuchintervall im rechten Bild.
        NCC_matrix = zeros(sp_L-2*s,l_d-r_d+1);

        % Iterieren durch alle Spalten von N_L bzw. N_R
        for X = 1:size(N_L,2)
            N_L_W = N_L(:,X);
            for Y = r_d:l_d
                if X+Y < 1 || X+Y > size(N_R,2)
                    % Pixel in Bild 2 liegt nicht innerhalb des Bildes mit
                    % mindestens Abstand s zum Bildrand
                    NCC_matrix(X,Y-r_d+1) = -1;
                else
                    % Korrelation mittels Skalarprodukt berechnen.
                    NCC_matrix(X,Y-r_d+1) = (1/(window_length^2-1))*(N_L_W'*N_R(:,X+Y));
                end
            end
        end
        [~, ind_max] = max(NCC_matrix,[],2);   

        % Spaltenvektor erstellen, der die PixelabstÃ¯Â¿Â½nde zwischen einem
        % Pixel im linken Bild und dem korrelierten Pixel im rechten Bild
        % enthÃ¯Â¿Â½lt. Erster Eintrag enthÃ¯Â¿Â½lt d vom Pixel mit der x-Koordinate
        % s+1 in der aktuellen Suchzeile im linken Bild.
        % d = m1 - m2 (Offset berÃ¼cksichtigen!)
        d_cor = off_d-(ind_max+r_d-1);
    end

function d_cor = SAD_line(line_number) % , min_cor
        % Function gibt nur d fÃ¯Â¿Â½r die Pixel des linkes Bildes mit mind.
        % Abstand s zum Bildrand zurÃ¯Â¿Â½ck => length(d_cor) = sp-2*s
               
        % Matrizen für stacked vectoren anlegen.
        N_L = zeros(window_length^2,sp_L-2*s);
        N_R = zeros(window_length^2,sp_R-2*s);
        % Normierte Fenster fÃ¼r linkes Bild berechnen:
        for k = s+1:sp_L-s
             % Quadratisches Fenster extrahieren.
             L_W = double(I_L_rec(line_number-s:line_number+s,k-s:k+s));
             % Mittelwert subtrahieren.
%              L_W = L_W - mean(L_W(:));
             % Teilen durch die Standardabweichung
%              L_W = L_W/std(L_W(:));
             % Normierte Matrix als Vektor speichern.
             idx = k-s;
             N_L(:,idx) = L_W(:);
        end
        % analog fÃ¼r rechtes Bild
        for k = s+1:sp_R-s
            % Suche im rechten Bild eine Zeile Ã¼ber der Zeile des linken
            % Bilds
             R_W = double(I_R_rec(line_number-s:...
                 line_number+s,k-s:k+s));
%              R_W = R_W - mean(R_W(:));
%              R_W = R_W/std(R_W(:));
             idx = k-s;
             N_R(:,idx) = R_W(:);
        end
        % Suche fÃ¯Â¿Â½r jedes Pixel im linken Bild passendes Pixel im rechten
        % Bild in einem Intervall von max_d bis min_d Pixel.
        % NCC-matrix: Zeile entspricht x-Pixelkoordinate im linken Bild,
        % Spalte entspricht dem Liniensuchintervall im rechten Bild.
        SAD_matrix = zeros(sp_L-2*s,l_d-r_d+1);

        % Iterieren durch alle Spalten von N_L bzw. N_R
        for X = 1:size(N_L,2)
            N_L_W = N_L(:,X);
            for Y = r_d:l_d
                if X+Y < 1 || X+Y > size(N_R,2)
                    % Pixel in Bild 2 liegt nicht innerhalb des Bildes mit
                    % mindestens Abstand s zum Bildrand
                    SAD_matrix(X,Y-r_d+1) = Inf;
                else
                    % Korrelation mittels Skalarprodukt berechnen.
                    SAD_matrix(X,Y-r_d+1) = sum(abs(N_L_W-N_R(:,X+Y)));
                end
            end
        end
        [~, ind_max] = min(SAD_matrix,[],2);   

        % Spaltenvektor erstellen, der die PixelabstÃ¯Â¿Â½nde zwischen einem
        % Pixel im linken Bild und dem korrelierten Pixel im rechten Bild
        % enthÃ¯Â¿Â½lt. Erster Eintrag enthÃ¯Â¿Â½lt d vom Pixel mit der x-Koordinate
        % s+1 in der aktuellen Suchzeile im linken Bild.
        % d = m1 - m2 (Offset berÃ¼cksichtigen!)
        d_cor = off_d-(ind_max+r_d-1);
    end
end
