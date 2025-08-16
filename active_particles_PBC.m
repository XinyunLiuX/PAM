function dydt = active_particles_PBC(t,y,p)
   
  N = p.N;

  dydt = zeros(4*N,1);
  xi = y(    1:N  ); 
  xi = position_apply_PBC(xi,p.box_length);                      %%%%%%% wrapped position <<< Xinyun >>>
  yi = y(  N+1:2*N); 
  yi = position_apply_PBC(yi,p.box_length);                      %%%%%%% wrapped position <<< Xinyun >>>
  ui = y(2*N+1:3*N);
  vi = y(3*N+1:4*N);

  % Friction term
 %Ri = (1/p.tau)*(1- ( ui.^2 + vi.^2 )./( 1 + p.epsilon*cos(p.omega*t) ) );   % <<< OLD  Value >>>
  Ri = (1/p.tau)*(1- ( ui.^2 + vi.^2 ).*( 1 + p.epsilon*cos(p.omega*t) ) );   % Xinyun Reverse

  % Distance between particles 
  [deltaxi,deltayi] = point_distance(xi,yi,N);
  deltaxi = distance_apply_PBC(deltaxi, p.box_length);           %%%%%%% wrapped distance <<< Xinyun >>>
  deltayi = distance_apply_PBC(deltayi, p.box_length);           %%%%%%% wrapped distance <<< Xinyun >>>
  
  rij = sqrt(deltaxi.^2 + deltayi.^2);
 
  dydt(    1:N  ) = ui;
  dydt(  N+1:2*N) = vi;
  dydt(2*N+1:3*N) = Ri.*ui - sum(p.alpha.*deltaxi./(rij.^(p.p+1)),2);
  dydt(3*N+1:4*N) = Ri.*vi - sum(p.alpha.*deltayi./(rij.^(p.p+1)),2);
  
end


function [distx,disty] = point_distance(xi,yi,N)


% Given the (x,y) coordinates of N points, this fuction computes the 
% distance from each point to all other points

% Inputs
% xi - vector with x coordinates of N points
% yi - vector with y coordinates of N points
% N  - number of points

% Outputs
% distx - matrix (N x N-1). Each row i contains in its columns all the
% distances in x to all other j points

% disty - matrix (N x N-1). Each row i contains in its columns all the
% distances in y to all other j points


% Duplicate columns
X = kron(xi,ones(1,N));
Y = kron(yi,ones(1,N));

distx = X'-X;
disty = Y'-Y;

% Remove diagonal elements
% x
utri = triu(distx);
utri = utri(:,2:end); %upper triangular

ltri = tril(distx);
ltri = ltri(:,1:end-1); %lower triangular

distx = utri+ ltri;

% y
utri = triu(disty);
utri = utri(:,2:end); %upper triangular

ltri = tril(disty);
ltri = ltri(:,1:end-1); %lower triangular

disty = utri+ ltri;

end
