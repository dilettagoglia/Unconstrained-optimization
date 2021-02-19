%How to run SQD%
cd C:\Users\Diletta\Documents\MATLAB
load('data_to_run_SQD.mat')
n=2

%prima provo con questo x%
x=zeros(n,1);

fcontour( f, [-5 , 1], 'LineColor', 'k', 'LineWidth', 1);

%intanto guardo subito qual è la soluzione ottimale%
xStar

[ x , state ] = SQD(Q, q, x)
fStar = 0.5 * xStar' * Q * xStar + q' * xStar;
fStar

%stampo il valore ottimale in questo modo per confrontarlo con il risultato
%che ho ottenuto e verificare che siano quasi uguali (spprossimazione
%corretta della soluzione ottimale)
fprintf( '\t%1.8e' , fStar )

%ora provo con questo x%
x=2*rand(n,1)
%chiamo di nuovo fcontour e [x,state]%


%ora proviamo con un nuovo intervallo
%i level set vengono riscalati (zoom in verso il valore ottimale%
fcontour( f, [-2.6 , -2.4], 'LineColor', 'k', 'LineWidth', 1);
%e con un nuovo valore di x (compreso nell'intervallo)%
[ x , state ] = SQD(Q, q, [-2.4 ; -2.5])

%ora stampo anche la colonna relativa all'errore f(x) - f* per ogni iterazione
%e il rate (rapporto con l'errore all'iterazione precedente) per capire l'efficienza dell'algoritmo
[ x , state ] = SQD(Q, q, [-2.4 ; -2.5], fStar)