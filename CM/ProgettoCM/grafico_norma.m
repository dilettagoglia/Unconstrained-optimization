% DISEGNARE IL GRAFICO 2D DELLA FUNZIONE IN COORDINATE POLARI
% questo perchè il valore della norma non dipende dal raggio
A = [ 1 2  ; 3 4 ;  5 6 ];
Q = A'*A;
theta = linspace(0, 2*pi); % theta va da zero a 2pigreco
x = [cos(theta) ; sin(theta)];

% f = - ( (x'*Q*x) / (x'*x) ); 
% ma dato che (x'*x)= cos^2 + sin^2 = 1
% disegno il grafico solo del numeratore



    % ----------------------------------------------- DISEGNARE LEVEL SET
    contour(f) % OK 
    % fcontour(f) non funziona perchè vuole input  function_handle
    % -------------------------------------------------------------------
    
plot(theta, f)
xlabel({'theta', '0 < x < 2\pi'})
ylabel('f')
% Il punto di minimo nel grafico 2D sarà in y = -(norm(A)^2)
% ------------------------------------------------------------

% DISEGNARE GRAFICO 3D CARTESIANO 
% definisco gli assi x1 e x2
x1=cos(theta);
x2=sin(theta);
plot3(x1, x2, f)
xlabel('x1')
ylabel('x2')
zlabel('f')

% NB: visto solo sul piano xy (ovvero dall'alto o dal basso),
% il grafico è la circonferenza unitaria !!
% ------------------------------------------------------------

% DISEGNARE IL GRAFICO 2D CARTESIANO 
% visualizzare circonferenza unitaria e vettore x
b=x/norm(x)
plot(x1, x2)
axis equal
xlabel('x1')
ylabel('x2')
hold on
plot(b) % NON FUNZIONA RIVEDERE

