function [image_corrected] = gain_offset_correction_cdf(image)
%GAIN_CORRECTION_CDF Summary of this function goes here
%   Beschreibung: Ausnutzen des dynamischen Intensitaetbereiches von 1% bis
%   99%
[bin_counts edges_intensity] = histcounts(image(:));
num_total = sum(bin_counts);
cdf = zeros(1,size(bin_counts,2));
cdf(1) = bin_counts(1)/num_total;
idx_low = 1;
idx_high = size(bin_counts,2)';
for i = 2:(length(bin_counts))
 cdf(i) = cdf(i-1)+bin_counts(i)/num_total;
 if cdf(i) < 0.01 %untere Schwelle der Intensitaet
        idx_low = i; 
 end
 if cdf(i) < 0.99 %obere Schwelle der Intensitaet
         idx_high = i;
 end
end

intensity_low = floor(edges_intensity(idx_low));
intensity_high = floor(edges_intensity(idx_high));
if(intensity_low <10 && intensity_high >250)
    %Robustheit gegen mehrmaliges Anwenden und schon sehr dynamische Bilder
    image_corrected = uint8(image);
else
%Ausnutzen des vollen Intensitaetsbereiches durch Streckung
image_corrected = ...
    uint8(double(image)-intensity_low)*(255/(intensity_high-intensity_low));
end



