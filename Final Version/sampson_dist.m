function sd = sampson_dist(F, x1_pixel, x2_pixel)
    % Diese Funktion berechnet die Sampson Distanz basierend auf der
    % Fundamentalmatrix F
    e3 = [0 -1 0; 1 0 0; 0 0 0];
    sd = dot(x2_pixel,(F*x1_pixel)).^2./(dot((e3*F*x1_pixel),(e3*F*x1_pixel))+dot((x2_pixel'*F*e3)',(x2_pixel'*F*e3)'));
end