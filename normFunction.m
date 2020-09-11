function normFunction = NormFunction()
normFunction = @norm2;

function [v, varargout] = norm2(x)

type basic_matrix.txt;
A = readmatrix('basic_matrix.txt');
[a, b] = size(A);


%% Esempio di risultato con questa matrice: 
%% norma ottenuta con il mio vettore x: 22.8433
%% norma di questa matrice: 23.0914

%% Potrebbe anche non essere così male, visto che la prima norma la calcoliamo tramite norm(A*x)/norm(x) 
%% => quindi eseguendo alcune operazioni che si portano dietro un certo errore


if isempty( x )  % informative call
   v = 0;
   if nargout > 1 
      % che vettore iniziale dare??
      varargout{ 1 } = ones(b,1);
   end
else
   Q = A'* A;
   v = (x'*Q*x)/(x'*x);  % f(x)
   % v = - norm(A*x)/norm(x);c
   if nargout > 1
      g = ((2*x)*(x'*Q*x))/((x'*x)^2) - ((2*Q*x)/(x'*x));  %gradiente
      varargout{ 1 } = g; 
   end
end
end 

end
