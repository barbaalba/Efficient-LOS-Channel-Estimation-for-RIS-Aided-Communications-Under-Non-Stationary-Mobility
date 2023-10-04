function [x_t,y_t] = randomwalk(NumUE,RandomLength,Xmax,Ymax,SP)
% Theis function geneate the location of the user within the room boundary
%%%%%%%%%%%%%%%%%%% INPUT %%%%%%%%%%%%%%%%%%%%%%
% NumUE: number of users in the room
% RandomLength: length of the simulation 
% Xmax,Ymax: They establish the size of the room 
% SP: Speed of the mobile user
%%%%%%%%%%%%%%%%%%% OUTPUT %%%%%%%%%%%%%%%%%%%%
% x_t,y_t; The location of the user with respect to time

x_t = zeros(NumUE,RandomLength);
y_t = zeros(NumUE,RandomLength);
for m=1:NumUE
  for n = 1:RandomLength 
      if x_t(m,n) < (Xmax-SP) && x_t(m,n)>(-Xmax+SP)
          A = SP*sign(randn); % sign value determines direction of user movement
          x_t(m,n+1) = x_t(m,n) + A;
      else
          A = -SP * sign(x_t(m,n));
          x_t(m,n+1) = x_t(m,n) + A; % if next move go beyond boundary move the oposite direction
      end
      if y_t(m,n) < Ymax && y_t(m,n) > -Ymax
          A = SP*sign(randn); % sign value determines direction of user movement
          y_t(m,n+1) = y_t(m,n) + A;
      else
          A = -SP*sign(y_t(m,n));
          y_t(m,n+1) = y_t(m,n) + A; % if next move go beyond boundary move the oposite direction
      end
  end
end
