function [T1, R1, T2, R2, U, V]=TR_aus_E(E)
    % Diese Funktion berechnet die moeglichen Werte fuer T und R
    % aus der Essentiellen Matrix
    % Rotationsmatrizen
    ROT_pPI2 = [0 -1 0; 1 0 0; 0 0 1];
    ROT_mPI2 = [0 1 0; -1 0 0; 0 0 1];
    
    % Singulärwertzerlegung von E
    [U,S,V] = svd(E);

    
    if(round(det(U))== -1)
        U = U*diag([1 1 -1]);
    end
    if(round(det(V))==-1)
        V = V*diag([1 1 -1]);
    end
    
    
    R1 = U*ROT_pPI2'*V';
    R2 = U*ROT_mPI2'*V';
    T1_mat = U*ROT_pPI2*S*U';
    T2_mat = U*ROT_mPI2*S*U';
    
    
    % Rekonstruktion des Translationsvektors
    T1(1,:) = T1_mat(3,2);
    T1(2,:) = T1_mat(1,3);
    T1(3,:) = T1_mat(2,1);
    T2(1,:) = T2_mat(3,2);
    T2(2,:) = T2_mat(1,3);
    T2(3,:) = T2_mat(2,1);
end