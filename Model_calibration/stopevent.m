function [value,isterminal,direction] = stopevent(tt,xx,Ap,Aorifice,Cd,rho,speed,S)

% This file is used to reject all the runs simulations that take longer
% than the set time. This runs are normally not useful and can be rejected
% without any issues.

if toc > 80 % time in seconds
    value = 0;
    isterminal = 1;     % stop the integration
    direction = 0;      % all events
else
    value = 1;
    isterminal = 0;     % don't stop
    direction = 0;      % all events
end