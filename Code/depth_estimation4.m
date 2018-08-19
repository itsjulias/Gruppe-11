function depth_map4 = depth_estimation4(I_L_rec,I_R_rec, off_d)
%DEPTH_ESTIMATION2 mittels NCC Liniensuche

% I_rec und I2_rec haben selbe Anzahl an Zeilen, aber nicht unbedingt selbe
% Anzahl an Spalten
[zl_L,sp_L] = size(I_L_rec);
sp_R = size(I_R_rec,2);

konst_dev = 0;

disparity_map = zeros(sp_L,zl_L);

        max_d = -200;
        min_d = 500;

window_length = 9;
% Erlaubte minimale Distanz zum Bildrand.
s = (window_length-1)/2;

% Liniensuche zeilenweise
            % Suche im rechten Bild eine Zeile über/unter der Zeile des linken
            % Bilds => konst_dev
            % Überarbeiten für konst_dev > 0!
for i = s+1-konst_dev:zl_L-s
% for i = 800:1000
    i/(zl_L-2*s+konst_dev)
%     (i-800)/200
    disparity_map(s+1:sp_L-s,i) = NCC_line(i);
end

disparity_map_unscaled = disparity_map;
disparity_map = disparity_map_unscaled-(off_d-min_d);
disparity_map = disparity_map ./ max(disparity_map(:));
figure
imshow(disparity_map')
% Überarbeiten!
depth_map4 = 1./disparity_map;
depth_map4 = depth_map4./10;
figure;
imshow(depth_map4');

%save('depth4_test_9_9_full_minus10')


    function d_cor = NCC_line(line_number) % , min_cor
        % Function gibt nur d fï¿½r die Pixel des linkes Bildes mit mind.
        % Abstand s zum Bildrand zurï¿½ck => length(d_cor) = sp-2*s
               
        % Matrizen für stacked vectoren anlegen.
        N_L = zeros(window_length^2,sp_L-2*s);
        N_R = zeros(window_length^2,sp_R-2*s);
        % Normierte Fenster für linkes Bild berechnen:
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
        % analog für rechtes Bild
        for k = s+1:sp_R-s
            % Suche im rechten Bild eine Zeile über der Zeile des linken
            % Bilds
             R_W = double(I_R_rec(line_number-s+konst_dev:...
                 line_number+s+konst_dev,k-s:k+s));
             R_W = R_W - mean(R_W(:));
             R_W = R_W/std(R_W(:));
             idx = k-s;
             N_R(:,idx) = R_W(:);
        end
        % Suche fï¿½r jedes Pixel im linken Bild passendes Pixel im rechten
        % Bild in einem Intervall von max_d bis min_d Pixel.
        % NCC-matrix: Zeile entspricht x-Pixelkoordinate im linken Bild,
        % Spalte entspricht dem Liniensuchintervall im rechten Bild.
        NCC_matrix = zeros(sp_L-2*s,min_d-max_d+1);

        % Iterieren durch alle Spalten von N_L bzw. N_R
        for X = 1:size(N_L,2)
            N_L_W = N_L(:,X);
            for Y = max_d:min_d
                if X+Y < 1 || X+Y > size(N_R,2)
                    % Pixel in Bild 2 liegt nicht innerhalb des Bildes mit
                    % mindestens Abstand s zum Bildrand
                    NCC_matrix(X,Y-max_d+1) = -1;
                else
                    % Korrelation mittels Skalarprodukt berechnen.
                    NCC_matrix(X,Y-max_d+1) = (1/(window_length^2-1))*(N_L_W'*N_R(:,X+Y));
                end
            end
        end
        [~, ind_max] = max(NCC_matrix,[],2);   

        % Spaltenvektor erstellen, der die Pixelabstï¿½nde zwischen einem
        % Pixel im linken Bild und dem korrelierten Pixel im rechten Bild
        % enthï¿½lt. Erster Eintrag enthï¿½lt d vom Pixel mit der x-Koordinate
        % s+1 in der aktuellen Suchzeile im linken Bild.
        % d = m1 - m2 (Offset berücksichtigen!)
        d_cor = off_d-(ind_max+max_d-1);
    end

end
