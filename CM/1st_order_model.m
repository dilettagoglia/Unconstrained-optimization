x = 0:0.1:3;
f=sin(x).*exp(x); % funzione %
d=diff(f); % derivata (= gradiente in una dimensione) %
L1=f(1)+d(1).*(x-1); % modello del primo ordine %
plot(f)
hold on
plot(d)
hold on
plot(L1)

% se mi avvicino ad un punto %
