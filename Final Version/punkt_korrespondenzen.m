function Korrespondenzen = punkt_korrespondenzen(I1,I2,Mpt1,Mpt2,varargin)
    % In dieser Funktion sollen die extrahierten Merkmalspunkte aus einer
    % Stereo-Aufnahme mittels NCC verglichen werden um Korrespondenzpunktpaare
    % zu ermitteln.
    
    %% Input parser
    p = inputParser;
    valid_window_length = @(x) (x>1) & (mod(x,2) == 1) & (mod(x,1)==0);
    valid_min_corr = @(x) isnumeric(x) & (x>0) & (x<1);
    valid_do_plot = @(x) islogical(x);
    addOptional(p,'window_length',25,valid_window_length);
    addOptional(p,'min_corr',0.95,valid_min_corr);
    addOptional(p,'do_plot',false,valid_do_plot);
    
    parse(p, varargin{:});
    
    Im1 = double(I1);
    Im2 = double(I2);
    window_length = p.Results.window_length;
    min_corr = p.Results.min_corr;
    do_plot = p.Results.do_plot;
    
    %% Merkmalsvorbereitung
    [length_y_im1, length_x_im1] = size(Im1);
    [length_y_im2, length_x_im2] = size(Im2);
    
    % Alle Merkmalspunkte, die zu nahe am Rand liegen werden nicht
    % berücksichtigt/entfernt
    Mpt1 = Mpt1(:,(Mpt1(1,:)+floor(window_length/2)<=length_x_im1) & (Mpt1(1,:)-floor(window_length/2)>=1)...
                  & (Mpt1(2,:)+floor(window_length/2)<=length_y_im1) & (Mpt1(2,:)-floor(window_length/2)>=1));
    Mpt2 = Mpt2(:,(Mpt2(1,:)+floor(window_length/2)<=length_x_im2) & (Mpt2(1,:)-floor(window_length/2)>=1)...
                  & (Mpt2(2,:)+floor(window_length/2)<=length_y_im2) & (Mpt2(2,:)-floor(window_length/2)>=1));
    
    %% Normierung
    Mat_feat_1 = zeros(window_length^2,length(Mpt1(1,:)));
    Mat_feat_2 = zeros(window_length^2,length(Mpt2(1,:)));
    window_dist = floor(window_length/2);
    for i=1:length(Mpt1(1,:))
        % Bildausschnitt extrahieren
        coordinate_pixel = Mpt1(:,i);
        section = Im1((coordinate_pixel(2)-window_dist):(coordinate_pixel(2)+window_dist),...
                            (coordinate_pixel(1)-window_dist):(coordinate_pixel(1)+window_dist));
        mean_section = mean(section(:));
        sigma_section = sqrt(1/(numel(section)-1)*trace((section(:)-mean_section)'*(section(:)-mean_section)));
        norm_section = 1/sigma_section*(section-mean_section);
        Mat_feat_1(:,i) = norm_section(:);
        
    end
    for i=1:length(Mpt2(1,:))
        % Bildausschnitt extrahieren
        coordinate_pixel = Mpt2(:,i);
        section = Im2((coordinate_pixel(2)-window_dist):(coordinate_pixel(2)+window_dist),...
                            (coordinate_pixel(1)-window_dist):(coordinate_pixel(1)+window_dist));
        mean_section = mean(section(:));
        sigma_section = sqrt(1/(numel(section)-1)*trace((section(:)-mean_section)'*(section(:)-mean_section)));
        norm_section = 1/sigma_section*(section-mean_section);
        Mat_feat_2(:,i) = norm_section(:);
        
    end
        
    %% NCC Brechnung
    [length_y_mat_feat1, length_x_mat_feat1] = size(Mat_feat_1);
    %[length_y_mat_feat2, length_x_mat_feat2] = size(Mat_feat_2)
    traceWV = Mat_feat_2'*Mat_feat_1;
    
    NCC_matrix = 1/(length_y_mat_feat1-1)*traceWV;
    NCC_matrix = NCC_matrix.*(double(1-(NCC_matrix<min_corr)));
    [val, sorted_index] = sort(NCC_matrix(:),'descend');
    
    %% Korrespondenz
    [size_NCC_y, size_NCC_x] = size(NCC_matrix);

    hKorrespondenzen = [];
    for i=1:length(sorted_index)
        index_y = mod(sorted_index(i),size_NCC_y);
        if(index_y == 0)
            index_y = size_NCC_y;
        end
        index_x = ceil(sorted_index(i)/size_NCC_y);
        if (sum(NCC_matrix(:,index_x))~= 0) 
            p1 = Mpt1(:,index_x);
            p2 = Mpt2(:,index_y);
            hKorrespondenzen = [hKorrespondenzen [p1; p2]];
            NCC_matrix(:,index_x) = zeros(size_NCC_y,1);
        end
    end
    Korrespondenzen = hKorrespondenzen;
    
    if(do_plot)
        figure;
        im1 = imshow(I1); hold all;
        im2 = imshow(I2);
        im1.AlphaData = 0.5;
        im2.AlphaData = 0.5;
        
        for i=1:length(Korrespondenzen(1,:))
            plot(Korrespondenzen(1,i), Korrespondenzen(2,i),'*r')
            plot(Korrespondenzen(3,i), Korrespondenzen(4,i),'*g')
            plot([Korrespondenzen(1,i), Korrespondenzen(3,i)], ...
                 [Korrespondenzen(2,i), Korrespondenzen(4,i)],'-b')
        end
        
    end
end