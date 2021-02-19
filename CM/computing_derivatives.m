cd C:\Users\Diletta\Documents\MATLAB
syms x y
%definisco variabili simboliche
%ora definisco funzione (drop)
f(x,y) = - (1 + cos(12*sqrt(x*x + y*y))) / (0.5*(x*x + y*y) + 2)

%calcolo derivata di f in x
diff(f,x)

%calcolo derivata seconda
diff(f,x,2)

%hessiana ??
diff( diff( f,x ), y)

%altro metodo: adigator