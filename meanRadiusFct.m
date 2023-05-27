function [position,isterminal,direction] = meanRadiusFct(t,y)
  position = norm(y(1:3)) - 6.25; % The value that we want to be zero
  isterminal = 1;  % Halt integration 
  direction = 0;   % The zero can be approached from either direction
end