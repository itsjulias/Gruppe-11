function new_rectified_img = new_rectified_view(I1_rec,depth_map,K,T_x)
%NEW_RECTIFIED_VIEW Gibt rektifizierte Zwischenansicht zurück
% -----------------------------------------------------------------------
% !!! T_x zwischen 0 (linke Ansicht) und -1 (rechte Ansicht) wählen !!!
% -----------------------------------------------------------------------

% P2 = R*P1+T mit R=I und P1=lamda1*x1
% lambda2*x2 = lambda1*x1+T => lambda2 ergibt sich durch Normierung der
% z-Komponente von P2 auf 1.
% Analoge Rechnung für Pixelkoordinaten:
% => lambda2*x2_pixel = lambda1*K*R*inv(K)*x1_pixel+K*T <= term_1
% P2 = inv(K)*K*lambda2*x2 = inv(K)*lambda2*x2_pixel
% => lambda2 ergibt sich durch Normierung der
% z-Komponente von P2 auf 1 => lambda2 = (3. Zeile von inv(K))*term_1 
% lambda1 ist durch depth_map gegeben, K*R*inv(K) = I und K*T = T_pixel
% => lambda2*x2_pixel = lambda1*x1_pixel+T_pixel <= term_1

% T_x soll zwischen 0 und -1 liegen
% Nur T_x nötig, da nur Verschiebung in x-Richtung (T = [T_x; 0; 0])
% T_pixel = K*T = [K(1,1)*T_x; 0; 0]

[L_pixel_x,L_pixel_y] = meshgrid(1:size(I1_rec,2),1:size(I1_rec,1));

P2_x_pixel = depth_map'.*L_pixel_x+T_x;
% T(2) = 0 und T(3) = 0, aus leterem folgt, dass P2_z = depth_map'
P2_y_pixel = depth_map'.*L_pixel_y;
% lambda2 ergibt sich durch Normierung der z-Komponente von P2 auf 1
% Skalarprodukt von 3. Zeile von inv(K) mit
% [P2_x_pixel P2_y_pixel P2_z_pixel]
K_inv = inv(K);
lambda2 = K_inv(3,1).*(P2_x_pixel+K(1,1)*T_x) + K_inv(3,2).*P2_y_pixel +...
          K_inv(3,3).*depth_map'; 
% Normierung mit lambda2
NEW_pixel_x = P2_x_pixel./lambda2;
NEW_pixel_y = P2_y_pixel./lambda2;
% Valide Pixel bestimmen
p_valid = isfinite(NEW_pixel_x) & isfinite(NEW_pixel_y);
NEW_pixel_x = NEW_pixel_x(p_valid);
NEW_pixel_y = NEW_pixel_y(p_valid);
v = double(I1_rec(p_valid));

F = scatteredInterpolant(NEW_pixel_x(:),NEW_pixel_y(:),v);
fprintf('F berechnet')
new_rectified_img = uint8(F({1:size(I1_rec,2),1:size(I1_rec,1)}))';

end

