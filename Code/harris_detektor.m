function merkmale = harris_detektor(input_image, varargin)
    % In dieser Funktion soll der Harris-Detektor implementiert werden, der
    % Merkmalspunkte aus dem Bild extrahiert
    
    %% Input parser
    p = inputParser;
    valid_segment_length = @(x) (x>=1) & (mod(x,2) == 1) & (mod(x,1)==0);
    valid_k = @(x) isnumeric(x) & (x>=0) & (x<=1);
    valid_tau = @(x) isnumeric(x) & (x>0);
    valid_do_plot = @(x) islogical(x);
    valid_min_dist = @(x) isnumeric(x) & x>=1;
    valid_tile_size = @(x) isnumeric(x);
    valid_N = @(x) isnumeric(x);
    addOptional(p,'segment_length',15,valid_segment_length);
    addOptional(p,'k',0.05,valid_k);
    addOptional(p,'tau',1e6,valid_tau);
    addOptional(p,'do_plot',false,valid_do_plot);
    addOptional(p,'min_dist',20,valid_min_dist);
    addOptional(p,'tile_size',200,valid_tile_size);
    addOptional(p,'N',5,valid_N);
    parse(p,varargin{:})
    
    if(isscalar(p.Results.tile_size))
       tile_size = [p.Results.tile_size p.Results.tile_size];
    else
       tile_size = p.Results.tile_size;
    end
    
    
    segment_length = p.Results.segment_length;
    k = p.Results.k;
    tau = p.Results.tau;
    do_plot = p.Results.do_plot;
    min_dist= p.Results.min_dist;
    N = p.Results.N;
    
    %% Vorbereitung zur Feature Detektion
    % Pruefe ob es sich um ein Grauwertbild handelt
    if(ndims(input_image) ~=2)
       error('Image format has to be NxMx1'); 
    end
    
    image = double(input_image);
    % Approximation des Bildgradienten
    [Ix,Iy] = sobel_xy(image);
    
    % Gewichtung
    sigma = floor(segment_length/2)/3;
    w = 1/(sqrt(2*pi)*sigma)*exp(-(-floor(segment_length/2):1:floor(segment_length/2)).^2/(2*sigma^2));
    % Harris Matrix G
    W = w'*w;
    G11 = conv2(Ix.^2,W,'same');
    G22 = conv2(Iy.^2,W,'same');
    G12 = conv2(Ix.*Iy,W,'same');
    
    %% Merkmalsextraktion ueber die Harrismessung
    [length_vertical, length_horizontal] =size(input_image);
    %H = zeros(length_vertical,length_horizontal);
    
    H = G11.*G22-(G12.^2)-k*(G11+G22).^2;
%     for i=1:length_vertical
%         for j=1:length_horizontal
%             G = [G11(i,j) G12(i,j); G12(i,j) G22(i,j)];
%             H(i,j) = det(G)-k*(G(1,1)+G(2,2))^2;
%         end
%     end

    corners = zeros(length_vertical,length_horizontal);
    border = ceil(segment_length/2);
    %h_corners = max(H-tau,0);
    h_corners = H.*(1-double(H<tau));
    corners(border:end-border, border:end-border) = h_corners(border:end-border, border:end-border);
    
    %% Merkmalsvorbereitung
    [length_vertical,length_horizontal] = size(corners);
    null_rand_corners = zeros(length_vertical+2*min_dist, length_horizontal+2*min_dist);
    null_rand_corners(min_dist+1:end-min_dist,min_dist+1:end-min_dist) = corners;
    corners = null_rand_corners;
    %transposed_corners = corners';
    vector_merkmale = corners(:);
    %[row]=find(vector_merkmale~=0);
    %[height, width] = size(corners);
    [val, idx] = sort(vector_merkmale,'descend');
    sorted_index = idx(val~=0);
    
    %% Akkumulatorfeld
    [height, width] = size(input_image);
    AKKA = zeros(ceil(height/tile_size(1)),ceil(width/tile_size(2)));
    
    
    %% Merkmalsbestimmung mit Mindestabstand und Maximalzahl pro Kachel

    [length_y, length_x] = size(corners);
    new_merkmale = [];
    circ = cake(min_dist);
    coordinate_pixel = [ceil(sorted_index'/length_y); mod(sorted_index',length_y)];
    coordinate_tile = [ceil((coordinate_pixel(1,:)-min_dist)/tile_size(1));...
        ceil((coordinate_pixel(2,:)-min_dist)/tile_size(2))];
    y1 = (coordinate_pixel(1,:)-min_dist);
    y2 = (coordinate_pixel(1,:)+min_dist);
    x1 = (coordinate_pixel(2,:)-min_dist);
    x2 = (coordinate_pixel(2,:)+min_dist);

    for i=1:length(sorted_index)
        % Abfrage ob Kachel schon voll ist 
        if(AKKA(coordinate_tile(2,i),coordinate_tile(1,i))<N)
           value = corners(coordinate_pixel(2,i),coordinate_pixel(1,i));
           if(value~=0)
               corners((x1(i):x2(i)),(y1(i):y2(i))) = corners((x1(i):x2(i)),(y1(i):y2(i))).*circ;
               corners(coordinate_pixel(2,i),coordinate_pixel(1,i)) = value;
               new_merkmale = [new_merkmale coordinate_pixel(:,i)];
               AKKA(coordinate_tile(2,i),coordinate_tile(1,i)) = ...
                 AKKA(coordinate_tile(2,i),coordinate_tile(1,i))+1;
           end   
        end
    end
    merkmale = new_merkmale-[min_dist; min_dist];
    
    %% Plot
    if(do_plot)
        figure
        size_image = size(input_image);
        imshow(input_image); hold all;
        plot(merkmale(1,:),merkmale(2,:),'*r');
    end
    
        
end