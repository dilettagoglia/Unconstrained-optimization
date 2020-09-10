function [ x , status ] =  SDG( f , varargin )

%function [ x , status ] = SDG( f , x , eps , MaxFeval , m1 , m2 , astart ,
%                               tau , sfgrd , mina )
%
% Apply the classical Steepest Descent algorithm for the minimization of
% the provided function f, which must have the following interface:
%
%   [ v , g ] = f( x )
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
% - g (real, [ n x 1 ] real vector): this also depends on x. if x == []
%   this is the standard starting point from which the algorithm should
%   start, otherwise it is the gradient of f() at x (or a subgradient if
%   f() is not differentiable at x, which it should not be if you are
%   applying the gradient method to it).
%
% The other [optional] input parameters are:
%
% - x (either [ n x 1 ] real vector or [], default []): starting point.
%   If x == [], the default starting point provided by f() is used.
%
% - eps (real scalar, optional, default value 1e-6): the accuracy in the 
%   stopping criterion: the algorithm is stopped when the norm of the
%   gradient is less than or equal to eps. If a negative value is provided,
%   this is used in a *relative* stopping criterion: the algorithm is
%   stopped when the norm of the gradient is less than or equal to
%   (- eps) * || norm of the first gradient ||.
%
% - MaxFeval (integer scalar, optional, default value 1000): the maximum
%   number of function evaluations (hence, iterations will be not more than
%   MaxFeval because at each iteration at least a function evaluation is
%   performed, possibly more due to the line search).
% 
% Parameters for line search:
%
% - m1 (real scalar, optional, default value 0.01): first parameter of the
%   Armijo-Wolfe-type line search (sufficient decrease). Has to be in (0,1)
%
% - m2 (real scalar, optional, default value 0.9): typically the second
%   parameter of the Armijo-Wolfe-type line search (strong curvature
%   condition). It should to be in (0,1); if not, it is taken to mean that
%   the simpler Backtracking line search should be used instead
%
% - astart (real scalar, optional, default value 1): starting value of
%   alpha in the line search (> 0) 
%
% - tau (real scalar, optional, default value 0.9): scaling parameter for
%   the line search. In the Armijo-Wolfe line search it is used in the
%   first phase: if the derivative is not positive, then the step is
%   divided by tau (which is < 1, hence it is increased). In the 
%   Backtracking line search, each time the step is multiplied by tau
%   (hence it is decreased).
%
% - sfgrd (real scalar, optional, default value 0.01): safeguard parameter
%   for the line search. to avoid numerical problems that can occur with
%   the quadratic interpolation if the derivative at one endpoint is too
%   large w.r.t. the one at the other (which leads to choosing a point
%   extremely near to the other endpoint), a *safeguarded* version of
%   interpolation is used whereby the new point is chosen in the interval
%   [ as * ( 1 + sfgrd ) , am * ( 1 - sfgrd ) ], being [ as , am ] the
%   current interval, whatever quadratic interpolation says. If you
%   experiemce problems with the line search taking too many iterations to
%   converge at "nasty" points, try to increase this
%
% - mina (real scalar, optional, default value 1e-16): if the algorithm
%   determines a stepsize value <= mina, this is taken as an indication
%   that something has gone wrong (the gradient is not a direction of
%   descent, so maybe the function is not differentiable) and computation
%   is stopped. It is legal to take mina = 0, thereby in fact skipping this
%   test.
%
% Output:
%
% - x ([ n x 1 ] real column vector): the best solution found so far.
%
% - status (string): a string describing the status of the algorithm at
%   termination 
%
%   = 'optimal': the algorithm terminated having proven that x is a(n
%     approximately) optimal solution, i.e., the norm of the gradient at x
%     is less than the required threshold
%
%   = 'stopped': the algorithm terminated having exhausted the maximum
%     number of iterations: x is the bast solution found so far, but not
%     necessarily the optimal one
%
%   = 'error': the algorithm found a numerical error that prevents it from
%     continuing optimization (see mina above)
%
%{
 =======================================
 Author: Antonio Frangioni
 Date: 13-11-19
 Version 1.11
 Copyright Antonio Frangioni
 =======================================
%}

Plotf = true;  % if f and the trajectory have to be plotted when n = 2

% reading and checking input- - - - - - - - - - - - - - - - - - - - - - - -
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

if ~ isa( f , 'function_handle' )
   error( 'f not a function' );
end

if isempty( varargin ) || isempty( varargin{ 1 } )
   [ fStar , x ] = f( [] );
else
   x = varargin{ 1 };
   if ~ isreal( x )
      error( 'x not a real vector' );
   end

   if size( x , 2 ) ~= 1
      error( 'x is not a (column) vector' );
   end
   
   fStar = f( [] );
end

n = size( x , 1 ); % size of vector x

if length( varargin ) > 1
   eps = varargin{ 2 };
   if ~ isreal( eps ) || ~ isscalar( eps )
      error( 'eps is not a real scalar' );
   end
else
   eps = 1e-4;
end

if length( varargin ) > 2
   MaxFeval = round( varargin{ 3 } );
   if ~ isscalar( MaxFeval )
      error( 'MaxFeval is not an integer scalar' );
   end
else
   MaxFeval = 1000; % max number of iterations
end

if length( varargin ) > 3
   m1 = varargin{ 4 };
   if ~ isscalar( m1 )
      error( 'm1 is not a real scalar' );
   end
   if m1 <= 0 || m1 >= 1
      error( 'm1 is not in (0 ,1)' );
   end       
else
   m1 = 0.01;                       % primo parametro della A-W line search
end

if length( varargin ) > 4
   m2 = varargin{ 5 };            
   if ~ isscalar( m1 )
      error( 'm2 is not a real scalar' );
   end
else
   m2 = 0.9;                      % secondo parametro della A-W line search
end

AWLS = ( m2 > 0 && m2 < 1 );      % controllo se parametro m2 soddisfa i valori, altrimenti faccio backtraking e non AW

if length( varargin ) > 5
   astart = varargin{ 6 };
   if ~ isscalar( astart )
      error( 'astart is not a real scalar' );
   end
   if astart < 0
      error( 'astart must be > 0' );
   end       
else
   astart = 0.1; % di default era 1, l'ho modificata e l'alg converge con molte meno iterazioni
end

if length( varargin ) > 6
   tau = varargin{ 7 };
   if ~ isscalar( tau )
      error( 'tau is not a real scalar' );
   end
   if tau <= 0 || tau >= 1
      error( 'tau is not in (0 ,1)' );
   end       
else
   tau = 0.9;
end

if length( varargin ) > 7
   sfgrd = varargin{ 8 };
   if ~ isscalar( sfgrd )
      error( 'sfgrd is not a real scalar' );
   end
   if sfgrd <= 0 || sfgrd >= 1
      error( 'sfgrd is not in (0, 1)' );
   end
else
   sfgrd = 0.01;
end

if length( varargin ) > 9
   mina = varargin{ 10 };
   if ~ isscalar( mina )
      error( 'mina is not a real scalar' );
   end
   if mina < 0
      error( 'mina is < 0' );
   end       
else
   mina = 1e-16; % check if steps are becoming too small (defaul value set to machine precision)
end

% "global" variables- - - - - - - - - - - - - - - - - - - - - - - - - - - -
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

lastx = zeros( n , 1 );  % (x^i) last point visited in the line search % iterazione xi
lastg = zeros( n , 1 );  % (g^i) gradient of lastx % direzione 
func_eval = 1;               % f() evaluations count ("common" with LSs)

% initializations - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

fprintf( 'Gradient method\n');
if fStar > - Inf
   fprintf( 'func_eval \t rel gap \t || g(x) || \t rate');
   prevv = Inf;
else
   fprintf( 'func_eval \t f(x) \t || g(x) ||');
end
fprintf( '\t ls feval \t a*');
fprintf( '\n\n' );

[ v , g ] = f( x );
ng = norm( g );
if eps < 0
   ng0 = - ng;  % norm of first subgradient: why is there a "-"? ;-)
else
   ng0 = 1;     % un-scaled stopping criterion
end

% main loop - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

while true
 
    % output statistics - - - - - - - - - - - - - - - - - - - - - - - - - -
    
    if fStar > - Inf
       fprintf( '%4d\t%1.4e\t%1.4e' , func_eval , ...
                ( v - fStar ) / max( [ abs( fStar ) 1 ] ) , ng ); %errore relativo
       if prevv < Inf
          fprintf( '\t%1.4e' , ( v - fStar ) / ( prevv - fStar ) ); %decremento dell'errore ad ogni iterazione
       else
          fprintf( '\t\t\t' );           
       end
       prevv = v;
    else
       fprintf( '%4d\t%1.8e\t\t%1.4e' , func_eval , v , ng );
    end
   
    % stopping criteria - - - - - - - - - - - - - - - - - - - - - - - - - -

    if ng <= eps * ng0
       status = 'optimal';
       break;
    end
    
    if func_eval > MaxFeval
       status = 'stopped';
       break;
    end
    
    % compute step size - - - - - - - - - - - - - - - - - - - - - - - - - -

    phip0 = - ng * ng;

    if AWLS
       [ a , v ] = ArmijoWolfeLS( v , phip0 , astart , m1 , m2 , tau );
    else
       [ a , v ] = BacktrackingLS( v , phip0 , astart , m1 , tau );
    end

    % output statistics - - - - - - - - - - - - - - - - - - - - - - - - - -

    fprintf( '\t\t%1.4e' , a );
    fprintf( '\n' );

    if a <= mina
       status = 'error';
       break;
    end

    % compute new point - - - - - - - - - - - - - - - - - - - - - - - - - -

    % possibly plot the trajectory
    if n == 2 && Plotf
       PXY = [ x ,  lastx ];
       line( 'XData' , PXY( 1 , : ) , 'YData' , PXY( 2 , : ) , ...
             'LineStyle' , '-' , 'LineWidth' , 2 ,  'Marker' , 'o' , ...
             'Color' , [ 0 0 0 ] );
       pause;
    end
    
    x = lastx;
  
    % update gradient - - - - - - - - - - - - - - - - - - - - - - - - - - -

    g = lastg;
    ng = norm( g );

    % iterate - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

end

% end of main loop- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% inner functions - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%calcola la restrizione della funzione phi lungo la direzione
function [ phi , phip ] = f2phi( alpha ) % gli passo alpha e mi restituisce sia funzione (phi) che derivata prima (phip) in alpha 
 % phi( alpha ) = f( x - alpha * g )
 % phi'( alpha ) = < \nabla f( x - alpha * g ) , g >

   lastx = x - alpha * g;
   [ phi , lastg ] = f( lastx );
   phip = - g' * lastg;
   func_eval = func_eval + 1;
end

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

function [ a , phia ] = ArmijoWolfeLS( phi0 , phip0 , as , m1 , m2 , tau )

% performs an Armijo-Wolfe Line Search.
%
% phi0 = phi( 0 ), phip0 = phi'( 0 ) < 0
%
% as > 0 is the first value to be tested: if phi'( as ) < 0 then as is
% divided by tau < 1 (hence it is increased) until this does not happen
% any longer
%
% m1 and m2 are the standard Armijo-Wolfe parameters; note that the strong
% Wolfe condition is used
%
% returns the optimal step and the optimal f-value


% FIRST PHASE -------------------------------------------------------

lsiter = 1;  % count iterations of first phase
while func_eval <= MaxFeval  % stopping condition 
   [ phia , phips ] = f2phi( as );

   if ( phia <= phi0 + m1 * as * phip0 ) && ... %Armijo
                                           ( abs( phips ) <= - m2 * phip0 ) %strong Wolfe
      fprintf( '\t%2d' , lsiter );
      a = as; %ho trovato lo step per cui valgono sia A che W e quindi mi fermo
      return;  % Armijo + strong Wolfe satisfied, we are done

   end
   if phips >= 0 %se derivata positiva mi fermo
      break;
   end
   as = as / tau; %altrimenti incremento lo step e riprovo
   lsiter = lsiter + 1;
end    

fprintf( '\t%2d ' , lsiter );
lsiter = 1;  % count iterations of second phase

% SECOND PHASE --------------------------------------------------

am = 0; %step iniziale dove la derivata è negativa
a = as; %step dove la derivata è positiva
phipm = phip0;

while ( func_eval <= MaxFeval ) && ( ( as - am ) ) > mina && ( phips > 1e-12 )

   % a = ( am * phips - as * phipm ) / ( phips - phipm ); % quadratic interpolation formula: gives me the point in the middle
   % compute the new value by safeguarded quadratic interpolation
   a = max( [ am + ( as - am ) * sfgrd min( [ as - ( as - am ) * sfgrd  a ] ) ] );

   % compute phi( a )
   [ phia , phip ] = f2phi( a );

   if ( phia <= phi0 + m1 * a * phip0 ) && ( abs( phip ) <= - m2 * phip0 )
      break;  % Armijo + strong Wolfe satisfied, we are done
   end

   % restrict the interval based on sign of the derivative in a
   if phip < 0
      am = a;
      phipm = phip;
   else
      as = a;
      if as <= mina
         break;
      end
      phips = phip;
   end
   lsiter = lsiter + 1;
end

fprintf( '%2d' , lsiter );

end

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

function [ as , phia ] = BacktrackingLS( phi0 , phip0 , as , m1 , tau )

% performs a Backtracking Line Search.
%
% phi0 = phi( 0 ), phip0 = phi'( 0 ) < 0
%
% as > 0 is the first value to be tested, which is decreased by
% multiplying it by tau < 1 until the Armijo condition with parameter
% m1 is satisfied
%
% returns the optimal step and the optimal f-value

lsiter = 1;  % count ls iterations
while func_eval <= MaxFeval && as > mina
   [ phia , ~ ] = f2phi( as );
   if phia <= phi0 + m1 * as * phip0  % Armijo satisfied
      break;                          % done
   end
   as = as * tau;
   lsiter = lsiter + 1;
end

fprintf( '\t%2d' , lsiter );

end

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

end  % the end- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -




