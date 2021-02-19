function [ x , status ] =  SDQ( Q , q , x , varargin )

%function [ x , status ] = SDQ( Q , q , x , fStar , eps , MaxIter )
%
% Apply the Steepest Descent algorithm with exact line search to the
% quadratic function
%
%   f(x) = 1/2 x^T Q x + q x
%
% Input:
%
% - Q ([ n x n ] real symmetric matrix, not necessarily positive
%   semidefinite): the Hessian (quadratic part) of f
%
% - q ([ n x 1 ] real column vector): the linear part of f
%
% - x ([ n x 1 ] real column vector): the point where to start the
%   algorithm from.
%
% - fStar (real scalar, optional, default value Inf): optimal value of f.
%   if a non-Inf value is provided it is used to print out stasistics about
%   the convergence speed of the algorithm
%
% - eps (real scalar, optional, default value 1e-6): the accuracy in the 
%   stopping criterion: the algorithm is stopped when the norm of the
%   gradient is less than or equal to eps
%
% - MaxIter (integer scalar, optional, default value 1000): the maximum
%   number of iterations
%
% Output:
%
% - x ([ n x 1 ] real column vector): either the best solution found so far
%   (possibly the optimal one) or a direction proving the problem is
%   unbounded below, depending on case
%
% - status (string): a string describing the status of the algorithm at
%   termination 
%
%   = 'optimal': the algorithm terminated having proven that x is a(n
%     approximately) optimal solution, i.e., the norm of the gradient at x
%     is less than the required threshold
%
%   = 'unbounded': the algorithm terminated having proven that the problem
%     is unbounded below: x contains a direction along which f is
%     decreasing to - Inf, either because f is linear along x and the
%     directional derivative is not zero, or because x is a direction with
%     negative curvature
%
%   = 'stopped': the algorithm terminated having exhausted the maximum
%     number of iterations: x is the best solution found so far, but not
%     necessarily the optimal one
%
%{
 =======================================
 Author: Antonio Frangioni
 Date: 08-10-17
 Version 0.10
 Copyright Antonio Frangioni
 =======================================
%}

Plotf = true;  % if f and the trajectory have to be plotted when n = 2

% reading and checking input- - - - - - - - - - - - - - - - - - - - - - - -
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

if ~ isreal( Q )
   error( 'Q not a real matrix' );
end

n = size( Q , 1 );

if n <= 1
   error( 'Q is too small' );
end

if n ~= size( Q , 2 )
   error( 'Q is not square' );
end

if ~ isreal( q )
   error( 'q not a real vector' );
end

if size( q , 2 ) ~= 1
   error( 'q is not a (column) vector' );
end

if size( q , 1 ) ~= n
   error( 'q size does not match with Q' );
end

if ~ isreal( x )
   error( 'x not a real vector' );
end

if size( x , 2 ) ~= 1
   error( 'x is not a (column) vector' );
end

if size( x , 1 ) ~= n
   error( 'x size does not match with Q' );
end

if isempty( varargin )
   fStar = - Inf;
   if ~ isreal( fStar ) || ~ isscalar( fStar )
      error( 'fStar is not a real scalar' );
   end
else
   fStar = varargin{ 1 };
end

if length( varargin ) > 1
   eps = varargin{ 2 };
   if ~ isreal( eps ) || ~ isscalar( eps )
      error( 'eps is not a real scalar' );
   end
   if eps < 0
      error( 'eps can not be negative' );
   end
else
   eps = 1e-6;
end

if length( varargin ) > 2
   MaxIter = round( varargin{ 3 } );
   if ~ isscalar( MaxIter )
      error( 'MaxIter is not an integer scalar' );
   end
else
   MaxIter = 1000;
end

% initializations - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

fprintf( 'Gradient method\n');
fprintf( 'iter\tf(x)\t\t\t||nabla f(x)||');
if fStar > - Inf
   fprintf( '\tf(x) - f*\trate');
   prevv = Inf;
end
fprintf( '\n\n' );

i = 1;

% main loop - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

while true
    % compute function value and gradient - - - - - - - - - - - - - - - - -

    v = 0.5 * x' * Q * x + q' * x;
    g = Q * x + q;
    ng = norm( g );

    % output statistics - - - - - - - - - - - - - - - - - - - - - - - - - - -
  
    fprintf( '%4d\t%1.8e\t\t%1.4e' , i , v , ng );
    if fStar > - Inf
       fprintf( '\t%1.4e' , v - fStar );
       if prevv < Inf
          fprintf( '\t%1.4e' , ( v - fStar ) / ( prevv - fStar ) );
       end
       prevv = v;
    end
    fprintf( '\n' );
   
    % stopping criteria - - - - - - - - - - - - - - - - - - - - - - - - - -
    if ng <= eps
       status = 'optimal';
       break;
    end
    
    if i > MaxIter
       status = 'stopped';
       break;
    end
    
    % compute step size - - - - - - - - - - - - - - - - - - - - - - - - - -
    % meanwhile, check if f is unbounded below

    den = g' * Q * g;
    
    if den <= 1e-12
       % this is actually two different cases:
       %
       % - g' * Q * g = 0, i.e., f is linear along g, and since the
       %   gradient is not zero, it is unbounded below
       %
       % - g' * Q * g < 0, i.e., g is a direction of negative curvature for
       %   f, which is then necessarily unbounded below
       %
       %fprintf( '%1.4e\n' , den );

       status = 'unbounded';
       break;
    end
        
    t = ng^2 / den;  % stapsize
    
    % compute new point - - - - - - - - - - - - - - - - - - - - - - - - - -

    % possibly plot the trajectory
    if n == 2 && Plotf
       PXY = [ x ,  x - t * g ];
       line( 'XData' , PXY( 1 , : ) , 'YData' , PXY( 2 , : ) , ...
             'LineStyle' , '-' , 'LineWidth' , 2 ,  'Marker' , 'o' , ...
             'Color' , [ 0 0 0 ] );
       pause;
    end
    
    x = x - t * g;
  
    % iterate - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    i = i + 1;

end

% end of main loop- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

end  % the end- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
