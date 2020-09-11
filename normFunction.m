function normFunction = NormFunction()
normFunction = @norm2;

function [v, varargout] = norm2(x)

type basic_matrix.txt;
A = readmatrix('basic_matrix.txt');
[a, b] = size(A);

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