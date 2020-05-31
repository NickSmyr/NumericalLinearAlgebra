%%%%%%%%%%%% ERGASIA 9 %%%%%%%%%%%%%%%%%%%
function PageRank
  
  A = randi([0 1] ,500 , 500);
  % Summing the columns
  columnSums = sum(A , 2);
  A = A ./ columnSums;
  
  % ? is currently row stochastic because the line sum to 1
  
  A = A';
  % Now A is column stochastic
  %plotGoogleMatrixEigenvalues(A);
  question1(A);
  question2(A);
  
  question3(A);
  question4();
end

function question3(A)
  %% Erwthma 3 %%
  disp("-------------Question 3---------------")
  disp("A = 0.85")
  alpha = 0.85;
  findTop10PagesPowerMethod(A , alpha);
  
  disp("A = 0.5")
  alpha = 0.5;
  findTop10PagesPowerMethod(A , alpha);
  
  disp("A = 0.25")
  alpha = 0.25;
  findTop10PagesPowerMethod(A , alpha);
  
  disp("A = 0")
  alpha = 0;
  findTop10PagesPowerMethod(A , alpha);
end
function question1(A)
  plotGoogleMatrixEigenvalues(A);
  hold off;
end
function question2(A)
  %% Erwthma 2 %%
  disp("-------------Question 2---------------")
  disp("A = 0.85")
  G = GoogleMatrix(A , 0.85);
  findTop10Pages(G);
  
  disp("A = 0.5")
  G = GoogleMatrix(A , 0.5);
  findTop10Pages(G);
  
  disp("A = 0.25")
  G = GoogleMatrix(A , 0.25);
  findTop10Pages(G);
  
  disp("A = 0")
  G = GoogleMatrix(A , 0.0);
  findTop10Pages(G);
end

function question4()
  load Purdue77587.mat
  % Loading the file fills variablre A with a sparse matrix
  % Because of the powermethod the error will be ((a * l2) / l1)^k 
  % where l1 is the largest eigenvalue
  % in magnitute and l2 the second largest eigenValue in magnitute
  
  % Therefore an approximate formula for the second largest eigenValue
  % is l2 = (10 ^ (- 6 / numIter) ) / a
  disp("-------------Question 4---------------")  
  disp("A = 0.85")
  alpha = 0.85;
  findTop10PagesPowerMethodWithUrls(A , alpha , url);
  %Niter = 66
  %approximate l2 = 0.95427
  
  disp("A = 0.5")
  alpha = 0.5;
  findTop10PagesPowerMethodWithUrls(A , alpha , url);
  %Niter = 17
  %approximate l2 =  0.88734
  disp("A = 0.25")
  alpha = 0.25;
  findTop10PagesPowerMethodWithUrls(A , alpha , url);
  %Niter = 9
  %approximate l2 =  0.86177
  disp("A = 0")
  alpha = 0;
  findTop10PagesPowerMethodWithUrls(A , alpha , url);
end

function findTop10Pages(G)
  %% Finds the top 10 pages for input google matrix G using the build in 
  %% eig procedure
  n = size(G  , 1);
  [ S , D] = eig(G);
  spectralRadiusIndex = -1;
  for i = 1:n
    if abs(D(i,i) - 1) < 10 ^ -10
        spectralRadiusIndex = i;
    end
  end
  eigenVector = S( : , spectralRadiusIndex);
  normalizedEigenVector = eigenVector / norm1(eigenVector);
  normalizedEigenVector = normalizedEigenVector * sign(normalizedEigenVector(1));

  sortedEigVector = sort(normalizedEigenVector , 'descend') ;
  disp("Top 10 elements in the eigenVector ")
  disp(sortedEigVector(1 : min(10  ,size(sortedEigVector)(1))));
  
end

function findTop10PagesPowerMethod(A , alpha)
  %% Finds the top 10 pages for input connectivity matrix A using the power method
  
  
  [ eigenVector , eigenValue ] = powerMethodForGoogleMatrix(A , alpha);
  
  normalizedEigenVector = eigenVector / norm1(eigenVector);
  normalizedEigenVector = normalizedEigenVector * sign(normalizedEigenVector(1));

  sortedEigVector = sort(normalizedEigenVector , 'descend') ;
  disp("Top 10 elements in the eigenVector ")
  disp(sortedEigVector(1 : min(10  ,size(sortedEigVector)(1))));
  
end

function findTop10PagesPowerMethodWithUrls(A , alpha , urls)
  %% Finds the top 10 pages for input connectivity matrix A using the power method
  
  
  [ eigenVector , eigenValue ] = powerMethodForGoogleMatrix(A , alpha);
  
  normalizedEigenVector = eigenVector / norm1(eigenVector);
  normalizedEigenVector = normalizedEigenVector * sign(normalizedEigenVector(1));
  %% Each element in the eigenvecto will not have the index of the original page
  vectorPlusUrls = [ normalizedEigenVector , (1:size(A)(1))' ];
  %% In order to sort in descending order in octave this must be done
  %% because direction is not implemented in sort Rows
  vectorPlusUrls(: , 1 ) = -1 * vectorPlusUrls(: , 1 );
  vectorPlusUrls =  sortrows(vectorPlusUrls);
  vectorPlusUrls(: , 1 ) = -1 * vectorPlusUrls(: , 1 );
  
  
  disp("Top 10 elements in the eigenVector ")
  for i = 1 : 10
    rank = vectorPlusUrls(i , 1);
    url = urls(vectorPlusUrls(i , 2));
    printf("rank : %f url: %s\n" , rank , cell2mat(url) );
  end
end

function plotGoogleMatrixEigenvalues(A)
  %% Erwthma 1 %%
  n = size(A)(1);
  m = size(A)(2);
  e = ones(n , 1);
  
  alpha = 0.85
  G = GoogleMatrix(A  , alpha);
  plotEigenValues(G , "rx" , "a = 0.85")
  plotCircle(alpha , "r--");
  
  alpha = 0.5
  G = GoogleMatrix(A  , alpha);
  plotEigenValues(G , "go" , "a = 0.5")
  plotCircle(alpha , "g--");
  
  alpha = 0.25
  G = GoogleMatrix(A  , alpha);
  plotEigenValues(G , "b^" , 'a = 0.25')
  plotCircle(alpha , "b--");
  
  alpha = 0;
  G = GoogleMatrix(A  , alpha);
  plotEigenValues(G , "k^" , 'a = 0')
  plotCircle(alpha , "k--");
  legend
  
end
function [x , l] = powerMethodForGoogleMatrix(A , alpha)
  n = size(A)(1);
  m = size(A)(2);
  
  e = ones(n , 1);
  
  x = rand( n , 1);
  x = x / norm1(x);
  error = 500;
  iterations = 0;
  while error >= 10^-6
    xnew = (alpha * (A * x))  + ( ((1-alpha) / n ) * e * (e' * x) );
    l = norm1(xnew);
    % The error is the difference in norm value of two subsequent vectors
    error = norm1(xnew - x);
    iterations += 1;
    x = xnew / norm1(xnew);
  end
  printf("Finished after %d iterations\n" , iterations);
end
function plotEigenValues(A , FMT , DisplayName)
  eigenvals = eig(A);
  real_z1 = real(eigenvals);
  imag_z1 = imag(eigenvals);
  
  hold on
  plot(real_z1,imag_z1 , FMT , 'DisplayName',DisplayName);
  hold off
  
end
function G = GoogleMatrix(A , alpha)
  %% Returns the google matrix for the input matrix A %%
  n = size(A)(1);
  e = ones(n , 1);
  G = (alpha * A)  + ( ((1-alpha) / n ) * e * e');
end

function plotCircle(radius ,FMT)
  hold on
  th = 0:pi/50:2*pi;
  xunit = radius * cos(th);
  yunit = radius * sin(th);
  h = plot(xunit, yunit , FMT);
  hold off
end

function b = norm1(x)
  % Returns the physical norm 1 of input vector x
  b = sum(abs(x));
end
