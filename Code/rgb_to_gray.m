function gray_image = rgb_to_gray(input_image)
    % Diese Funktion soll ein RGB-Bild in ein Graustufenbild umwandeln. Falls
    % das Bild bereits in Graustufen vorliegt, soll es direkt zurueckgegeben werden.
    if(ndims(input_image) == 3)
        img = double(input_image);
        gray_image = uint8(img(:,:,1)*0.299 + img(:,:,2)*0.587 + img(:,:,3)*0.114);
    else
        gray_image = input_image;
    end
end