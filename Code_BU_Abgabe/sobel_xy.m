function [Fx, Fy] = sobel_xy(input_image)
    % In dieser Funktion soll das Sobel-Filter implementiert werden, welches
    % ein Graustufenbild einliest und den Bildgradienten in x- sowie in
    % y-Richtung zurueckgibt.
    sobel_matrix_x = [1 0 -1; 2 0 -2; 1 0 -1]; %horizontal
    sobel_matrix_y = sobel_matrix_x'; %vertikal
    
    Fx = conv2(input_image,sobel_matrix_x,'same');
    Fy = conv2(input_image,sobel_matrix_y,'same');

end