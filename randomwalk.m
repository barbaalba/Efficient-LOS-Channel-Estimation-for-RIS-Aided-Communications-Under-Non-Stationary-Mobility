function [x_t,y_t] = randomwalk(NumUE,RandomLength)
%N = 1000; % Length of the x-axis, also known as the length of the random walks.
%M = 2; % The amount of random walks.
x_t = zeros(NumUE,RandomLength);
y_t = zeros(NumUE,RandomLength);
for m=1:NumUE
  for n = 1:RandomLength % Looping all values of N into x_t(n).
    A = sign(randn); % Generates either +1/-1 depending on the SIGN of RAND.
    x_t(m,n+1) = x_t(m,n) + A;
    A = sign(randn); % Generates either +1/-1 depending on the SIGN of RAND.
    y_t(m,n+1) = y_t(m,n) + A;
  end
end
