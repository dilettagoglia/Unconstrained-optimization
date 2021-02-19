

A=randn(3,3); %matrice 3x3 di valori random%
[V, D]=eig(A); %eigenvalue decomposition%
               %D � quella matrice diagonale che chiamiamo lambda%
v = V(:,1); %prendo la prima colonna di V e la memorizzo in v%
A*v, v*D(1,1); %controllo che questi due risultati siano uguali%

%%%%%%%%%
B = 1/3 * A; %ha tutti autovalori < 1 (in modulo)%

C=eye(3); %identit�%
for k = 1:20
C = C*B
end %posso notare che C diventa sempre pi� piccola, 
    %infatti il limite per k-->infinito di C^k � zero
    %(se tutti gli autovalori sono < 1 in modulo)