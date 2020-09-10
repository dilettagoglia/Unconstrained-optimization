% last update: 10-09-2020 (Diletta)

function F = norma()

% Produces a cell array of a function handler, useful to run the unconstrained
% optimization algorithm of estimating the 2-norm for a possibly rectangular matrix A.
%
% The function in the array has the following interface:
%
%   [ v , varargout ] = f( x )
%
% Input:
%
% - x is either a [ n x 1 ] real (column) vector denoting the input of
%   f(), or [] (empty).
%
% Output:
%
% - v (real, scalar): if x == [] this is the best known lower bound on
%   the unconstrained global optimum of f(); it can be -Inf if either f()
%   is not bounded below, or no such information is available. If x ~= []
%   then v = f(x).
%
% - g (real, [ n x 1 ] real vector) is the first optional argument. This
%   also depends on x. if x == [] this is the standard starting point of an
%   optimization algorithm, otherwise it is the gradient of f() at x, or a
%   subgradient if f() is not differentiable at x.
%
% - H (real, [ n x n ] real matrix) is the first optional argument. This
%   must only be specified if x ~= [], and it is the Hessian of f() at x.
%   If no such information is available, the function throws error.
%{
 =======================================
 Author: Diletta Goglia, Lisa Lavorati
 Date: 10-09-2020
 Version 1.0
 Copyright D. Goglia, L. Lavorati
 =======================================
%}
    F = {[1]}; % cell( 1 , 1 );
    F{1}=@Norma;
    
%--------------------------------------------------------------------------

    function [ v , varargout ] = Norma( x ) % [ f(x) , nabla f(x) ]

    % a=3; % il valore dovrebbe essere random: a=randi(10); 
    b=2; 
    %A=rand(a,b) % matrice rettangolare (m x n) di valori random
     
    A = [ 1 2  ; 3 4 ;  5 6 ];

    if isempty( x )  % informative call
       v = 0;   

       if nargout > 1 % nargout returns the number of function output arguments specified in the call to the currently executing function
          varargout{ 1 } = ones(b, 1); %varargout  is an output variable in a function definition statement that enables the function to return any number of output arguments
       end
    else

        Q = A'*A; % Q è quadrata e simmetrica
        % x = sym('x', [n 1]) % x target vector (n x 1) --> x(1), x(2), ... , xn sono le variabli della funzione
        v = - ( (x'*Q*x) / (x'*x) );  % f(x)
        % v = - norm(A*x)/norm(x);

       if nargout > 1
         %calcolo derivate parziali in ogni variabile x e le memorizzo nel gradiente
         g = zeros( 2 , 1 );
         % g = ( ( 2*x*(x'*Q*x) ) / (x'*x)^2 ) - ( (2*Q*x) / (x'*x) );
         g( 1 ) = (2*x(1)*(x(1)*(35*conj(x(1)) + 44*conj(x(2))) + x(2)*(44*conj(x(1)) + 56*conj(x(2)))))/(x(1)*conj(x(1)) + x(2)*conj(x(2)))^2 - (70*x(1) + 88*x(2))/(x(1)*conj(x(1)) + x(2)*conj(x(2)));
         g( 2 ) = (2*x(2)*(x(1)*(35*conj(x(1)) + 44*conj(x(2))) + x(2)*(44*conj(x(1)) + 56*conj(x(2)))))/(x(1)*conj(x(1)) + x(2)*conj(x(2)))^2 - (88*x(1) + 112*x(2))/(x(1)*conj(x(1)) + x(2)*conj(x(2)));
         % [m,n] = size(g); %controllo le dimensioni del gradiente (1 x n)

         varargout{ 1 } = g;  % gradiente di  f(x)
         
         if nargout > 2
             % syms x1 x2
             % x = [x1;x2]
             % H = hessian(v, [x1,x2]);
             H = zeros( 2 , 2 );
             H( 1 , 1 ) = (2*(x1*(35*conj(x1) + 44*conj(x2)) + x2*(44*conj(x1) + 56*conj(x2))))/(x1*conj(x1) + x2*conj(x2))^2 - 70/(x1*conj(x1) + x2*conj(x2)) - (2*(x1 + conj(x1))^2*(x1*(35*conj(x1) + 44*conj(x2)) + x2*(44*conj(x1) + 56*conj(x2))))/(x1*conj(x1) + x2*conj(x2))^3 + (2*(x1 + conj(x1))*(35*x1 + 44*x2 + 35*conj(x1) + 44*conj(x2)))/(x1*conj(x1) + x2*conj(x2))^2;
             H( 1 , 2 ) = ((x2 + conj(x2))*(35*x1 + 44*x2 + 35*conj(x1) + 44*conj(x2)))/(x1*conj(x1) + x2*conj(x2))^2 - 88/(x1*conj(x1) + x2*conj(x2)) + ((x1 + conj(x1))*(44*x1 + 56*x2 + 44*conj(x1) + 56*conj(x2)))/(x1*conj(x1) + x2*conj(x2))^2 - (2*(x1 + conj(x1))*(x2 + conj(x2))*(x1*(35*conj(x1) + 44*conj(x2)) + x2*(44*conj(x1) + 56*conj(x2))))/(x1*conj(x1) + x2*conj(x2))^3;
             H( 2 , 2 ) = (2*(x1*(35*conj(x1) + 44*conj(x2)) + x2*(44*conj(x1) + 56*conj(x2))))/(x1*conj(x1) + x2*conj(x2))^2 - 112/(x1*conj(x1) + x2*conj(x2)) - (2*(x2 + conj(x2))^2*(x1*(35*conj(x1) + 44*conj(x2)) + x2*(44*conj(x1) + 56*conj(x2))))/(x1*conj(x1) + x2*conj(x2))^3 + (2*(x2 + conj(x2))*(44*x1 + 56*x2 + 44*conj(x1) + 56*conj(x2)))/(x1*conj(x1) + x2*conj(x2))^2;
             H( 2 , 1 ) = H( 1 , 2 );
             varargout{ 2 } = H;       % hessiana di f(x)
     
         end

       end
    end

    end 
%--------------------------------------------------------------------------
end 
