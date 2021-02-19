function TF = TestFunctions()

%function TF = TestFunctions()
%
% Produces a cell array of function handlers, useful to test unconstrained
% optimization algorithms.
%
% Each function in the array has the following interface:
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
%
% The current list of functions is the following:
%
%  1 Standard 2x2 PSD quadratic function with nicely conditioned Hessian.
%
%  2 Standard 2x2 PSD quadratic function with less nicely conditioned
%    Hessian.
%
%  3 Standard 2x2 PSD quadratic function with Hessian having one zero
%    eigenvalue.
%
%  4 Standard 2x2 quadratic function with indefinite Hessian (one positive
%    and one negative eigenvalue)
%
%  5 Standard 2x2 quadratic function with "very elongated" Hessian (a 
%    very small positive minimum eigenvalue, the other much larger)
%
%  6 the 2-dim Rosenbrock function
%
%  7 the "six-hump camel" function
%
%  8 the Ackley function
%
%  9 a 2-dim nondifferentiable function coming from Lasso regularization
%
%{
 =======================================
 Author: Antonio Frangioni
 Date: 08-11-18
 Version 1.01
 Copyright Antonio Frangioni
 =======================================
%}

TF = cell( 1 , 1 );
TF{ 1 } = @rosenbrock;

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

function [ v , varargout ] = rosenbrock( x )
% rosenbrock's valley-shaped function
% syms x y
% f = @(x, y) 100 * ( y - x^2 )^2 + ( x - 1 )^2
%
% diff( f , x )
% 2 * x - 400 * x * ( - x^2 + y ) - 2
%
% diff( f , y )
% - 200 * x^2 + 200 * y
%
% diff( f , x , 2 )
% 1200 * x^2 - 400 * y + 2
%
% diff( f , y , 2 )
% 200
%
% diff( f , x , y )
% -400 * x

if isempty( x )  % informative call
   v = 0;
   if nargout > 1
      varargout{ 1 } = [ -1 ; 1 ];
   end
else
   v = 100 * ( x( 2 ) - x( 1 )^2 )^2 + ( x( 1 ) - 1 )^2;  % f(x)
   if nargout > 1
      g = zeros( 2 , 1 );
      g( 1 ) = 2 * x( 1 ) - 400* x( 1 ) * ( x( 2 ) - x( 1 )^2 ) - 2;
      g( 2 ) = - 200 * x( 1 )^2 + 200 * x( 2 );
      
      varargout{ 1 } = g;  % \nabla f(x)
      if nargout > 2
         H = zeros( 2 , 2 );
         H( 1 , 1 ) = 1200 * x( 1 )^2 - 400 * x( 2 ) + 2;
         H( 2 , 2 ) = 200;
         H( 2 , 1 ) = -400 * x( 1 );
         H( 1 , 2 ) = H( 2 , 1 );
         varargout{ 2 } = H;       % \nabla^2 f(x)
      end
   end
end
end  % rosenbrock


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

end