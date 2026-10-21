function d_y = getGeneralizedOutputVector(sys,x_u,n_x,t)
% function d_y = getGeneralizedOutputVector(sys,x_u,n_x,t)
%
% Evaluates the generalized output vector (state derivatives and outputs)
% of the Simulink model sys at the working point x_u at time t.
%
% Inputs:
% sys       Name of the Simulink model (without file extension) as string
% x_u		Generalized input vector: column vector of states x followed by
%           inputs u
% n_x       Number of states. This is needed to disassemble x_u into x and u
% t         Time in seconds
%
% Outputs:
% d_y		Generalized output vector: column vector of state derivatives d
%           followed by outputs y

% Calculate outputs and derivatives at a working point (x,u).
y = feval (sys, t, x_u(1:n_x), x_u(n_x+1:end), 'outputs');
d = feval (sys, t, x_u(1:n_x), x_u(n_x+1:end), 'derivs');

d_y = [d; y];

end