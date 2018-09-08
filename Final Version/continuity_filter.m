function filtered_map = continuity_filter(map,thres,window_length)
%CONTINUITY_FILTER Function betrachtet die benachbarten Pixel innerhalb
%eines Fenster mit Seitenlänge window_length. Weicht der Mittelwert der
%benachbarten Pixel mehr als der Schwellwert 'thres' ab, wird das aktuelle
%Pixel als Ausreißer detektiert und durch den Mittelwert ersetzt.
%Andernfalls bleibt der Pixelwert unverändert.
filtered_map = map;

        % Abstand Pixel zu Fensterrand (Fensterlänge soll ungerade sein)
        s = (window_length-1)/2;
        
count = 0;
for m = s+1:size(map,1)-s
    for n = s+1:size(map,2)-s
        win = map(m-s:m+s,n-s:n+s);
        mean_rahmen = (sum(win(:))-win(s+1,s+1))/(window_length^2-1);
        mean_dif = mean_rahmen - win(s+1,s+1);
        if mean_dif > thres
            filtered_map(m,n) = mean_rahmen;
            count = count +1;
        end
    end
end
% count

