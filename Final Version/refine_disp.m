function disp_map = refine_disp(I1_rec,I2_rec,hom_map,interv_search_left_matrix,interv_search_right_matrix,off_d,window_length,L_R_base)
%REFINE_DISP Refine the disparity map

% I1_rec und I2_rec haben selbe Anzahl an Zeilen, aber nicht unbedingt selbe
% Anzahl an Spalten
[zl,spL] = size(I1_rec);
spR = size(I2_rec,2);

% Erlaubte minimale Distanz zum Bildrand.
s = (window_length-1)/2;

% I1_rec,I2_rec und hom_map sollen die selbe Größe haben, damit später ein
% linearer Index definiert werden kann, der für beide Bilder gültig ist.
if spR > spL
    I1_rec = [I1_rec,inf(zl,spR-spL)];
    hom_map = [hom_map,zeros(zl,spR-spL)];
    interv_search_left_matrix = [interv_search_left_matrix,nan(zl,spR-spL)];
    interv_search_right_matrix = [interv_search_right_matrix,nan(zl,spR-spL)];
elseif spR < spL
    I2_rec = [I2_rec,inf(zl,spL-spR)];
end
sp = max(spL,spR);

% Später transponieren
disp_map = nan(sp,zl);

% Sicherstellen dass hom_map keine Werte > 0 enthält die zu nah am Rand
% liegen.
% Linken Rand füllen
hom_map(:,1:s) = 0;
% Rechten Rand füllen
hom_map(:,end-s+1:end) = 0;
% Oberen Rand füllen
hom_map(1:s,:) = 0;
% Unteren Rand füllen
hom_map(end-s+1:end,:) = 0;

%% Korrespondenzsuche mittels SAD
% Alle Fenster im zweiten Bild als stacked vector in der dritten Dimension
% speichern. Transponieren um später Liniensuche mit Hilfe der linearen
% Indizes durchzuführen.
%% Fenster als stacked vector ermitteln
L_W_cell = cell(sp,zl);
R_W_cell = cell(sp,zl);
for k = s+1:zl-s
    for l = s+1:sp-s
        % Linkes Bildfenster nur berechnen falls keine homogene Fläche
        % vorliegt.
        if hom_map(k,l) > 0
            if interv_search_left_matrix(k,l) >= interv_search_right_matrix(k,l)
                % Linkes Bildfenster nur berechnen, falls Disparität nicht
                % bereits eindeutig.
                hom_map(k,l) = 0;
            else
                L_W = double(I1_rec(k-s:k+s,l-s:l+s));
                L_W_cell{l,k} = L_W(:);
            end
        end
        R_W = double(I2_rec(k-s:k+s,l-s:l+s));
        R_W_cell{l,k} = R_W(:);
    end
end
%% Ränder füllen
% Fülle leere Ränder von R_W_cell
len_sv = window_length^2;
% Oberer Rand
for k = 1:s
    for l = 1:sp
        R_W_cell{l,k} = inf(len_sv,1);
    end
end
% Unteren Rand
for k = zl-s+1:zl
    for l = 1:sp
        R_W_cell{l,k} = inf(len_sv,1);
    end
end
% Linken Rand
for k = s+1:zl-s
    for l = 1:s
        R_W_cell{l,k} = inf(len_sv,1);
    end
end
% Rechten Rand
for k = s+1:zl-s
    for l = sp-s+1:sp
        R_W_cell{l,k} = inf(len_sv,1);
    end
end

%% Spaltensuche statt Zeilensuche
% Indizes der Stellen finden, die nicht als homogene Fläche detektiert
% wurden. (Transponieren damit Zeilensuche zu Spaltensuche mit Hile des
% linearen Index wird)
idx_hom = find(hom_map');

% Fensterausschnitte sind nun spaltenweise gespeichert. (Ignoriere
% Bildränder, da Liniensuchbereich bereits durch niedrigere Auflösung
% vorgeben ist.)
interv_search_left_matrix = interv_search_left_matrix';
interv_search_right_matrix = interv_search_right_matrix';

%% Bildgrenzen einhalten
% Für Pixel oben links oder unten rechts im Bild kann es sein, dass der
% Suchbereich außerhalb des Bildes liegt. Diese Fälle sollen zunächst extra
% behandelt werden.
k = 1;
start_idx = idx_hom(k)+interv_search_left_matrix(idx_hom(k));
while start_idx < 1
    idx = idx_hom(k);
    % Schleife kann höchstens für ein paar Pixel der obersten Zeile des
    % Bildes erreicht werden. In diesem Fall nehme Mittelwert der
    % bisherigen min./max. Dispariäten.
    idx_sp = round((interv_search_right_matrix(idx)+interv_search_left_matrix(idx))/2);
    if strcmp(L_R_base,'R_base')
        disp_map(idx) = off_d+idx_sp;
    else
        disp_map(idx) = off_d-idx_sp;
    end
    k = k+1;
    start_idx = idx_hom(k)+interv_search_left_matrix(idx_hom(k));
end

m = length(idx_hom);
end_idx = idx_hom(m)+interv_search_right_matrix(idx_hom(m));
max_idx = zl*sp;
while end_idx > max_idx
    idx = idx_hom(k);
    % Schleife kann höchstens für ein paar Pixel der obersten Zeile des
    % Bildes erreicht werden. In diesem Fall nehme Mittelwert der
    % bisherigen min./max. Dispariäten.
    idx_sp = round((interv_search_right_matrix(idx)+interv_search_left_matrix(idx))/2);
    if strcmp(L_R_base,'R_base')
        disp_map(idx) = off_d+idx_sp;
    else
        disp_map(idx) = off_d-idx_sp;
    end
    m = m-1;
    end_idx = idx_hom(m)+interv_search_right_matrix(idx_hom(m));
end

%% Eigentlicher SAD-Algorithmus
if strcmp(L_R_base,'R_base')
    for i = k:m
        idx = idx_hom(i);
        interv_left = interv_search_left_matrix(idx);
        interv_right = interv_search_right_matrix(idx);
        start_idx = idx+interv_left;
        end_idx = idx+interv_right;
        % Indizes des Suchbereichs bestimmen
        search_interv = start_idx:end_idx;
        [~,idx_minSAD] = min(sum(abs(repmat(L_W_cell{idx},1,length(search_interv))-...
            [R_W_cell{search_interv}])));
        disp_map(idx) = off_d+interv_left+idx_minSAD-1;
    end
else
    for i = k:m
        idx = idx_hom(i);
        interv_left = interv_search_left_matrix(idx);
        interv_right = interv_search_right_matrix(idx);
        start_idx = idx+interv_left;
        end_idx = idx+interv_right;
        % Indizes des Suchbereichs bestimmen
        search_interv = start_idx:end_idx;
        [~,idx_minSAD] = min(sum(abs(repmat(L_W_cell{idx},1,length(search_interv))-...
            [R_W_cell{search_interv}])));
        disp_map(idx) = off_d-(idx_minSAD+interv_left-1);
    end
end

disp_map = disp_map';

if spR > spL
    disp_map = disp_map(:,1:spL);
end

end

