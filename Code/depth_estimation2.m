function depth_map2 = depth_estimation2(I1_rec,I2_rec)
%DEPTH_ESTIMATION2 mittel NCC Liniensuche

[zl sp] = size(I1_rec);

disparity_map = zeros(sp,zl);

window_length = 9;
% Erlaubte minimale Distanz zum Bildrand.
s = (window_length-1)/2;

% Liniensuche zeilenweise
% for i = s+1:zl-s
for i = 700:1000
    (i-700)/300
    disparity_map(s+1:sp-s,i) = NCC_line(i);
end
% �berarbeiten!
disparity_map_unscaled = disparity_map;
disparity_map = disparity_map -550;
% disparity_map = disparity_map;
disparity_map = disparity_map ./ max(disparity_map(:));
figure
imshow(disparity_map')
% �berarbeiten!
depth_map2 = 1./disparity_map;
depth_map2 = depth_map2./10;
figure;
imshow(depth_map2');

save('depth2_test_9_9_V3')


    function d_cor = NCC_line(line_number) % , min_cor
        % Function gibt nur d f�r die Pixel des linkes Bildes mit mind.
        % Abstand s zum Bildrand zur�ck => length(d_cor) = sp-2*s
        
        % Suchintervall [-500,200] evt. noch anpassen!
        
        % Matrizen f�r stacked vectoren anlegen.
        N_L = zeros(window_length^2,sp-2*s);
        N_R = zeros(window_length^2,sp-2*s);
        for k = s+1:sp-s
             % Quadratisches Fenster extrahieren.
             L_W = double(I1_rec(line_number-s:line_number+s,k-s:k+s));
             R_W = double(I2_rec(line_number-s:line_number+s,k-s:k+s));
             % Mittelwert subtrahieren.
             L_W = L_W - mean(L_W(:));
             R_W = R_W - mean(R_W(:));
             % Teilen durch die Standardabweichung
             L_W = L_W/std(L_W(:));
             R_W = R_W/std(R_W(:));
             % Normierte Matrix als Vektor speichern.
             idx = k-s;
             N_L(:,idx) = L_W(:);
             N_R(:,idx) = R_W(:);
        end
        % Suche f�r jedes Pixel im linken Bild passendes Pixel im rechten
        % Bild in einem Intervall von -550 bis +200 Pixel.
        % (V1 mit -600 bis 300)
        % NCC-matrix: Zeile entspricht x-Pixelkoordinate im linken Bild,
        % Spalte entspricht dem Liniensuchintervall im rechten Bild.
        
        max_d = -1350;
        min_d = -550;
        
        NCC_matrix = zeros(sp-2*s,min_d-max_d+1);
        

        % Iterieren durch alle Spalten von N_L bzw. N_R
        for X = 1:size(N_L,2)
            for Y = max_d:min_d
                if X+Y < 1 %|| X+Y > size(N_L,2)
                    % Pixel in Bild 2 liegt nicht innerhalb des Bildes mit
                    % mindestens Abstand s zum Bildrand
                    NCC_matrix(X,Y-max_d+1) = -1;% NCC_matrix(X,Y+1351) = -1;
                else
                    % Korrelation mittels Skalarprodukt berechnen.
                    % NCC_matrix(X,Y+1351) = (1/(window_length^2-1))*(N_R(:,X)'*N_L(:,X+Y));
                    NCC_matrix(X,Y-max_d+1) = (1/(window_length^2-1))*(N_L(:,X)'*N_R(:,X+Y));
                end
            end
        end
        [~, ind_max] = max(NCC_matrix,[],2);   
%         [V, ind_max] = max(NCC_matrix,[],2);
%         % Checke, ob min_cor eingehalten wird.
%         for m = 1:numel(V)
%             if V(m) < min_cor
%                 if m == 1
%                     % 1. Pixel im Zeilensuchbereich, linkes Pixel daneben
%                     % enth�lt null, daher suche n�chstes Pixel in rechte
%                     % Richtung bei dem V(m) > min_cor und �bernehme d.
%                     ind_max(m) = ind_max(find(V>min_cor,1));
%                 else
%                     % Falls nicht 1. Pixel, �bernehme d vom n�chsten
%                     % linksglegenen Pixel
%                     ind_max(m) = ind_max(m-1);
%                 end
%             end
%         end
        % Spaltenvektor erstellen, der die Pixelabst�nde zwischen einem
        % Pixel im linken Bild und dem korrelierten Pixel im rechten Bild
        % enth�lt. Erster Eintrag enth�lt d vom Pixel mit der x-Koordinate
        % s+1 in der aktuellen Suchzeile im linken Bild.
        % d = m1 - m2;
        d_cor = -(ind_max+max_d-1);
    end

end

