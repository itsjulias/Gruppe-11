function [features_img_1,features_img_2] = extract_1D_features(img1,img2,varargin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%% Input parser
p = inputParser;
valid_segment_length = @(x) (x>=1) & (mod(x,2) == 1) & (mod(x,1)==0);
valid_tau = @(x) isnumeric(x) & (x>0);
valid_do_plot = @(x) islogical(x);
valid_min_dist = @(x) isnumeric(x) & x>=1;
valid_N = @(x) isnumeric(x) & x>0;
valid_filter_kernel = @(x) (mod(x,1) == 0) & (x>0);

addOptional(p,'segment_length',15,valid_segment_length);
addOptional(p,'tau',0.5,valid_tau);
addOptional(p,'do_plot',false,valid_do_plot);
addOptional(p,'min_dist',20,valid_min_dist);
addOptional(p,'N',20,valid_N);
addOptional(p,'size_filter_kernel',10,valid_filter_kernel);
parse(p,varargin{:})



segment_length = p.Results.segment_length;
size_filter_kernel = p.Results.size_filter_kernel;
tau = p.Results.tau;
do_plot = p.Results.do_plot;
min_dist= p.Results.min_dist;
N = p.Results.N;

% Bias-Gain-Anpassung der Bilder, damit in etwa dieselben Merkmale gefunden
% werden
[img1_corrected] = gain_offset_correction_cdf(img1);
[img2_corrected] = gain_offset_correction_cdf(img2);

% Konvertierung der Bilder in double-Werte
img1_corrected = double(img1_corrected);
img2_corrected = double(img2_corrected);

% Gewichtung
sigma = floor(size_filter_kernel/2)/3;
w = 1/(sqrt(2*pi)*sigma)*exp(-(-floor(size_filter_kernel/2):1:floor(size_filter_kernel/2)).^2/(2*sigma^2));

% Für jedes Bild (img1 & img2) wird Intensität in x-Richtung untersucht.
% Dabei sollen Kanten entdeckt werden.

features_img_1 = [];
features_img_2 = [];

% Step 1: Gradientenberechnung mit anschließender Filterung durch
% Gaußfenster
[dI1_x, dI1_y] = sobel_xy(img1_corrected,'optimal',true);
[dI2_x, dI2_y] = sobel_xy(img2_corrected,'optimal',true);

% Entferne Bildränder, da diese nicht im Gradienten erscheinen sollen
dI1_x = dI1_x(2:end-1,2:end-1);
dI1_y = dI1_y(2:end-1,2:end-1);
dI2_x = dI2_x(2:end-1,2:end-1);
dI2_y = dI2_y(2:end-1,2:end-1);


%Step 1 (cont'd): Filterung mit Gaußfenster
dI1_filt = conv2(dI1_x.^2+dI1_y.^2,w'*w,'same');
dI2_filt = conv2(dI2_x.^2+dI2_y.^2,w'*w,'same');

for i=1:size(dI1_filt,1)
    
    %Step 2: Randmerkmale werden entfernt. Jede Zeile wird auf maximalen
    %Gradientenwert normiert, sodass Erkennung von Kanten durch Faktor in
    %Intervall [0, 1] möglich ist 
    dI1_filt(i,1:ceil(segment_length/2)) = zeros(1,ceil(segment_length/2));
    dI2_filt(i,1:ceil(segment_length/2)) = zeros(1,ceil(segment_length/2));
    dI1_filt(i,(end-floor(segment_length/2)):end) = ...
        zeros(1,ceil(segment_length/2));
    dI2_filt(i,(end-floor(segment_length/2)):end) = ...
        zeros(1,ceil(segment_length/2));

    dI1_filt(i,:) = dI1_filt(i,:);%./max(dI1_filt(i,:));
    dI2_filt(i,:) = dI2_filt(i,:);%./max(dI2_filt(i,:));

    %Step 3: Merkmale nach absteigender Reihenfolge sortieren und danach
    %Extrahieren der maximal zulässigen Merkmale pro Zeile (absteigend nach
    %Gradientenwert)  
    
    %Sortieren
    [values_features_1, idx_1] = sort(dI1_filt(i,:),'descend');
    idx_1 = idx_1(values_features_1>tau);
    values_features_1 = values_features_1(values_features_1>tau);
    
    [values_features_2, idx_2] = sort(dI2_filt(i,:),'descend');
    idx_2 = idx_2(values_features_2>tau);
    values_features_2 = values_features_2(values_features_2>tau);
    
    % minimaler Abstand soll garantiert werden
    matrix = zeros(1,length(idx_1));
    matrix(1) = 1;
    idx_valid_1 = [];
    while ~isempty(idx_1)
        idx_valid_1 = [idx_valid_1 idx_1(1)];
        matrix = (idx_1>(idx_1(1)+min_dist) | idx_1<(idx_1(1)-min_dist));
        idx_1 = idx_1(matrix);
    end
    
    matrix = zeros(1,length(idx_1));
    matrix(1) = 1;
    idx_valid_2 = [];
    while ~isempty(idx_2)
        idx_valid_2 = [idx_valid_2 idx_2(1)];
        matrix = (idx_2>(idx_2(1)+min_dist) | idx_2<(idx_2(1)-min_dist));
        idx_2 = idx_2(matrix);
    end
    
    %Clipping, wenn es zu viele Merkmale pro Zeile gibt
    if(length(idx_valid_1)>N)
        idx_valid_1 = idx_valid_1(1:N);
    elseif(length(idx_valid_2)>N)
        idx_valid_2 = idx_valid_1(1:N);
    end
        
    features_img_1 = [features_img_1,[idx_valid_1;i*ones(1,length(idx_valid_1))]];
    features_img_2 = [features_img_2,[idx_valid_2;i*ones(1,length(idx_valid_2))]];

end



if(do_plot)
    figure 
    imshow(uint8(img1_corrected));
    hold all
    plot(features_img_1(1,:),features_img_1(2,:),'x')
    
    figure
    imshow(uint8(img2_corrected));
    hold all
    plot(features_img_2(1,:),features_img_2(2,:),'x')
end

end

