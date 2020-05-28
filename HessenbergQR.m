function HessenbergQR
  [Q , R ] = qr(rand(10));
  randomVector = 60 * rand(10 , 1) - 30;
  D = diag(randomVector);
  A = Q * D * Q';
  
  disp("Executing builtin eig algorithm")
  tic
  eig(A);
  toc
  
  disp("Executing qr algorithm with initial hessenberg transform")
  tic
  upperTriang = HessenbergQRAlgorithm(A);
  toc
  
  disp("Executing qr algorithm with no preconditioning")
  tic
  upperTriang = baseQRAlgorithm(A);
  toc
  
  %%%% Apotelesmata:
  %%% Autes oi methodoi den katafernoun na sgklinoun gia 100 * 100 arxiko pinaka
  %%% akoma kai meta apo 25000 iterations
  %%% Wstoso gia mikroterous pinakes o algorithmos qr o opoios 
  %%% kanei ton pinaka hessenberg stin arxh einai pio grhgoros
  
  
  
end
function [FinalUpperTriangular] = baseQRAlgorithm(A)
  R = A;
  iterations = 0;
  while ~stopCriteria(R)
    [Q , R] = householderQR(R);
    R = R*Q;
    iterations +=1;
    if mod(iterations , 100) == 0 
      printf("Progress report iterations %d \n" , iterations);
      printf("Average abs %f \n" , averageAbsUnderDiagonal(R))
      printf("Max abs %f \n\n" , maxAbsUnderDiagonal(R))
    end
  
  end
  printf("Total iterations %d \n" , iterations);    
  FinalUpperTriangular = R;
end
function [FinalUpperTriangular] = HessenbergQRAlgorithm(A)
  % Performs the qr algorithm by first tranforming A to a hessenberg matrix
  % then performing the algorithm iterations with qr with givens
  [Qinit , T] = makeSimilarHessenberg(A);
  R = T;
  iterations = 0;
  while ~stopCriteria(R)
    [Q , R] = HessenbergGivensQR(R);
    R = R*Q;
    iterations +=1;
    if mod(iterations , 100) == 0 
      printf("Progress report iterations %d \n" , iterations);
      printf("Average abs %f \n" , averageAbsUnderDiagonal(R))
      printf("Max abs %f \n\n" , maxAbsUnderDiagonal(R))
    end
  end
  printf("Total iterations %d \n" , iterations);
  FinalUpperTriangular = R;
end


function [Q  , R] = HessenbergGivensQR(A)
  % Performs the qr decomposition of the input hessian matrix
  % With givens rotations so that  Q * R = A
  n = size(A)(1);
  m = size(A)(2);
  Q = eye(n);
  R = A;
  for i = 1 : m -1
    % Making the (i+1 . i ) element zero
    oldR = R;
    [ G , R ]  = givens(i , i+1 , R);
    %We need to multiply from the right because Q is the inverse of the product
    % of the givens rotations
    Q = Q * G';
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


function b = stopCriteria(A)
  % Returns true if every element under the main diagonal is smaller than
  % 10 ^ -5
  n = size(A)(1);
  m = size(A)(2);
  for j = 1:m
    for i = j+1: n
      if abs(A(i , j)) > 10^-5
        b = 0;
        return
      end
    end
  end
  b = 1;      
end
function b = averageAbsUnderDiagonal(A)
  n = size(A)(1);
  m = size(A)(2);
  b = 0;
  count = 0;
  for j = 1:m
    for i = j+1: n
      b += abs(A(i , j));
      count += 1;
    end
  end
  b = b / count;
end

function b = maxAbsUnderDiagonal(A)
  n = size(A)(1);
  m = size(A)(2);
  b = -999999;
  for j = 1:m
    for i = j+1: n
      if abs(A(i , j)) > b
        b = abs(A(i , j));
      end
    end
  end
end

function [Q  , T] = makeSimilarHessenberg(A)
  % Performs a similarity transform on matrix A 
  % to create hessenberg matrix T so that T = Q' * ? * Q
  
  % for every row (but not the last)
  n = size(A)(1);
  m = size(A)(2);
  Q = eye(n);
  T = A;
  for i = 1 : m -1 
    y = T(i + 1:end , i);
    ysize = n - i;
    c = norm(y);
    ut = y + sign(y(1)) * c * eye(ysize , 1);
    u = ut / norm(ut);    
    P = eye(ysize) - 2 * u * u';    
    P1 = blkdiag(eye(n - ysize) , P);
    Q =  Q * P1';
    T = P1 * T * P1';
  end
end

function [G , newA] = givens( k , j , A)
  % Performs a givens rotation on the matrix A so that the element
  % A(j , k ) is zero in the transformed matrix
  % Returns : G the givens rotation that gives newA when Applied to A
  % That Means G * A = newA
  s = - A(j , k) /( sqrt(A(k,k) ^ 2  + A(j,k)^2) );
  c =  A(k , k) / ( sqrt(A(k,k) ^ 2  + A(j,k)^2) );
  G = sparse(eye(size(A)(1)));
  G(k,k) = c;
  G(k,j) = -s;
  G(j,k) = s;
  G(j,j) = c;
  newA = A;
  % Only elements in the k , j rows are altered in the original Matrix A

  newA(k , : ) = c * A(k , :) - s * A(j , :);
  newA(j , : ) = s * A(k , :) + c * A(j , :); 
  
 
end

