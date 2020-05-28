function PageRank
  
  A = randi([0 1] ,10 , 10);
  % Summing the columns
  columnSums = sum(A , 2);
  A = A ./ columnSums;
  A
  
  plotGoogleMatrixEigenvalues(A);
  plotCircle(0.85 , "--r")
  plotCircle(0.5 , "--g")
end

function plotGoogleMatrixEigenvalues(A)
  A = A'
  alpha = 0
  n = size(A)(1);
  m = size(A)(2);
  nmmnsq = n * m - n^2
  e = ones(n , 1);

  G = GoogleMatrix(A'  , alpha);
  G
  [V , S] = eig(G);
  eigenvals = diag(S);
  
  max_eigenValue = max(eigenvals)
  real_z1 = real(eigenvals);
  imag_z1 = imag(eigenvals);
  
  plot(real_z1,imag_z1 , "or");
 
  real_z1 = real(eigenvals);
  imag_z1 = imag(eigenvals);
  
  plot(real_z1,imag_z1 , "or");
  
end

function G = GoogleMatrix(A , alpha)
  n = size(A)(1);
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
