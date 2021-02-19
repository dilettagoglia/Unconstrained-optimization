
% Copia e incolla i seguenti comandi:
A = [ 1 2  ; 3 4 ;  5 6 ];
Q = A'*A;
cd C:\Users\Diletta\Documents\MATLAB\Progetto % path to directory
F = norma();
[ x , state ] = SDG(F{1})

% DISEGNARE LEVEL SET --------------------------------
x1 = x(1)
x2 = x(2)
f = @(x1, x2) - ( ([x1 x2]*Q*[x1;x2])/([x1 x2]*[x1;x2]) ); % f = - ( (x'*Q*x) / (x'*x) );
fcontour( f, 'LineColor', 'k', 'LineWidth', 1)
% -----------------------------------------------------


% altri punti di partenza
r = rand(n , 1)
[ x , state ] = SDG(F{1}, r) 

% trovare vettori interessanti da cui partire 



% PROVARE CORRETTEZZA RISULTATO --------------------------------

%           Salvare vettore x ottenuto e calcolare:
%           le seguenti forme devono restutuire tutte valori uguali a norm(A)
%          - norm( A * x ) / norm( x ) 
%          - sqrt(x'*Q*x) / sqrt(x'*x)
%          - il primo valore di Sigma  -->  [U, Sigma, V]= svd(A);          


%       restituiti vettori diversi a seconda dei punti diversi da cui parto
%       MA posso verificare che sono tutti proporzionali normalizzandoli, ovvero
%       calcolando x/norm(x) per ogni vettore x restituito
%       e verificare che il risultato sia sempre lo stesso e sia esattamente
%       la prima colonna di V in -->  [U, Sigma, V]= svd(A);
%       (e che abbia norma 1)
%       dalla normalizzazione dei vettori deduco che 
%       il valore della norma non dipende dal raggio
% --------------------------------------------------------------


