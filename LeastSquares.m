function LeastSquares
  A = createMatrix
  b = [1 : 10]';
  disp("Condition of A")
  cond(A)
  x1 = normalEquations(A  , b);
  
  disp("Normal equations error")
  norm(A*x1 - b)
  
  x2 = LeastSquaresSVD(A  , b);
  
  disp("Svd error")
  norm(A*x2 - b)
  
  x3 = LeastSquaresHouseHolderQR(A , b);
  
  disp("Householder QR error")
  norm(A*x3 - b)
  
  x4 = LeastSquaresBaseGramSchmidtQR(A , b);
  
  disp("Base Gram Schmidt QR error")
  norm(A*x4 - b)
  
  x5 = LeastSquaresModifiedGramSchmidtQR(A , b);
  
  disp("Modified Gram Schmidt QR error")
  norm(A*x5 - b)
end

function x =  normalEquations(A , b)
  x = (A' * A)^-1 *   (A' * b);
end 
function x = LeastSquaresSVD(A , b)
  n = size(A)(1);
  m = size(A)(2);
  
  [U , S , V]  = svd(A);
  U = U(1 :n , 1:m);
  S = S(1 : m , 1 : m);
  x = V * S^-1 * U' * b;  
end 

function x = LeastSquaresHouseHolderQR(A , b)
  n = size(A)(1);
  m = size(A)(2);
  
  [Q , R] = householderQR(A);
  Q = Q(1: n , 1:m);
  R = R(1: m , 1:m);
  x = R^-1 * Q' * b;
  
end

function x = LeastSquaresBaseGramSchmidtQR(A , b)
  
  
  [Q , R] = basicGramSchmidt(A);
  x = R^-1 * Q' * b;
  
end

function x = LeastSquaresModifiedGramSchmidtQR(A , b)
  
  
  [Q , R] = basicGramSchmidt(A);
  x = R^-1 * Q' * b;
  
end

function [A]  = createMatrix()
  A = zeros(10 , 8);
  for i = 1 : 10
    for j = 1:8
      A(i , j) = 1 / ( i + j + 10 );
    end
  end

end

function [Q ,R] = householderQR(A)
  % Computes the qr decomposition with householder transformations
  n = size(A)(1);
  m = size(A)(2);
  Q = eye(n);
  R = A;
  
  for j = 1 : m
    y = A( j : end , j );
    ysize = size(y)(1);
    c = norm(y);
    ut = y + sign(y(1)) * c * eye(ysize , 1);
    u = ut / norm(ut);    
    P = eye(ysize) - 2 * u * u';
    P1 = blkdiag( eye(n - ysize) , P);
    R = P1 * R;
    Q  = Q * P1';
  end
end

function [Q , R] =  modifiedGramSchmidt(A)
  n = size(A)(1);
  m = size(A)(2);
  Q = zeros(n , m);
  R = zeros(m , m);
  
  for i = 1: m
    Q(1 : end , i) = A(1 : end , i);
    
    for j = 1 : i-1
      R(j , i)= A(1 : end , i)' * Q(1 : end , j);
      Q(1 : end , i) = Q(1 : end , i) - R(j , i) * Q(1 : end , j);
    end
    R(i,i) = norm(Q(1 : end , i));
    Q(1 : end , i) = Q(1 : end , i) / R(i,i); 
  end
end

function [Q , R] =  basicGramSchmidt(A)
  n = size(A)(1);
  m = size(A)(2);
  Q = zeros(n , m);
  R = zeros(m , m);
  
  for i = 1: m
    Q(1 : end , i) = A(1 : end , i);
    
    for j = 1 : i-1
      R(j , i)= Q(1 : end , i)' * Q(1 : end , j);
      Q(1 : end , i) = Q(1 : end , i) - R(j , i) * Q(1 : end , j);
    end
    R(i,i) = norm(Q(1 : end , i));
    Q(1 : end , i) = Q(1 : end , i) / R(i,i); 
  end
end

