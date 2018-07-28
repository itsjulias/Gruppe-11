function Cake = cake(min_dist)
    % Die Funktion cake erstellt eine "Kuchenmatrix", die eine kreisfoermige
    % Anordnung von Nullen beinhaltet und den Rest der Matrix mit Einsen
    % auffuellt. Damit koennen, ausgehend vom staerksten Merkmal, andere Punkte
    % unterdrueckt werden, die den Mindestabstand hierzu nicht einhalten. 
    idx_row = 0:min_dist;
    idx_col = 0:min_dist;
    B = (((idx_row)/min_dist)'.^2+((idx_col)/min_dist).^2)>1;
    Cake = [fliplr([fliplr(B(2:end,2:end)) B(2:end,:)]')';[fliplr(B(:,2:end)) B]];
    
    
end