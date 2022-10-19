function [x_t,y_t] = randomwalk(NumUE,RandomLength,Xmax,Ymax)
%N = 1000; % Length of the x-axis, also known as the length of the random walks.
%M = 2; % The amount of random walks.
x_t = zeros(NumUE,RandomLength);
y_t = zeros(NumUE,RandomLength);
for m=1:NumUE
  for n = 1:RandomLength % Looping all values of N into x_t(n).
      if x_t(m,n) < Xmax && x_t(m,n)>-Xmax
          A = sign(randn); % Generates either +1/-1 depending on the SIGN of RAND.
          x_t(m,n+1) = x_t(m,n) + A;
      else
          A = -sign(x_t(m,n));
          x_t(m,n+1) = x_t(m,n) + A;
      end
      if y_t(m,n) < Ymax && y_t(m,n) > -Ymax
          A = sign(randn); % Generates either +1/-1 depending on the SIGN of RAND.
          y_t(m,n+1) = y_t(m,n) + A;
      else
          A = -sign(y_t(m,n));
          y_t(m,n+1) = y_t(m,n) + A;
      end
  end
end
