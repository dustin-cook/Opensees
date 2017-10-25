%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function computes the time history response of an SDOF, for a period and a
% damping ratio.
%
% This was modified from the "elastic_Sa" function created originally by
% J.W. Baker.
%
% This automatically decreases the time step if T/10 < dT
% For hand notes on the verification of this processor, please see hand
% notes dated 8-16-05 in the ATC-63 folder.
%
% File was created by Jack Baker
% Send to Curt Haselton on 8-15-05
% File modified by Curt Haselton on 8-16-05
% File again modified by Curt Haselton to output response TH on 3-16-09
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[psuedoAccelerationTH, displacementTH, accelerationTH] = fn_sdof_th(T, z, ag, dt)

% INPUT 
% T     = period                            ( scalar, seconds )
% z     = damping ratio                     ( ", percent of critical damping )
% ag    = ground acceleration time history  ( length(ag) x 1 )
% Dtg   = time step of ag                   ( scalar )
% Dt    = analyis time step                 ( scalar )

% m = 1 (mass)

% OUTPUT 
%displacementTH = u(:,1);
%velocityTH = u(:,2);
%accelerationTH = u(:,3);


% Calculate displacement for the specified number of timesteps
% using the above information The equation of motion solved 
% is mu'' + cu' + ku = -m*ag(t), where m = 1, c = 2*xi*wn,
% k = wn^2, and p(t) = -u''_g, where u''_g is given

w = 2*pi ./ T;
k = w.^2;
gamma = 1/2;
beta = 1/4;
c = z .* (2*w);

k_hat = k + gamma/(beta*dt)*c + 1/(beta*dt^2);

a = 1/(beta*dt) + gamma/beta*c;
b = 1/(2*beta) + dt*(gamma/(2*beta) - 1)*c;

% if T/10 < dt, interpolate ground motion to refine analysis
if(T/10 < dt)
    tg = [ 0:(length(ag)-1) ]' * dt;
    new_dt = T/10;
    t = [ 0:new_dt:tg(end) ]';
    ag = interp1( tg, ag, t );
end


% displacement array with cols representing u, u', u''
u = zeros(length(ag),3);
u(1,3) = -ag(1);

for i = 1:(length(ag)-1)
    
   delta_p = -(ag(i+1) - ag(i)) + a*u(i,2) + b*u(i,3);
   delta_u = delta_p/k_hat;
   delta_u_dot = gamma/(beta*dt)*delta_u - gamma/beta*u(i,2) + dt*(1-gamma/(2*beta))*u(i,3);
   delta_u_dotdot = 1/(beta*dt^2)*delta_u - 1/(beta*dt)*u(i,2) - 1/(2*beta)*u(i,3);
   
   u(i+1,1) = u(i,1) + delta_u;
   u(i+1,2) = u(i,2) + delta_u_dot;
   u(i+1,3) = u(i,3) + delta_u_dotdot;
  
end

Sa = max(abs(u(:,1))) * (2*pi/T)^2;

displacementTH = u(:,1);
velocityTH = u(:,2);
accelerationTH = u(:,3);

psuedoAccelerationTH = u(:,1) .* (2*pi/T)^2;


