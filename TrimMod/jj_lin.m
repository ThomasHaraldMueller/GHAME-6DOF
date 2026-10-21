function jaco = jj_lin (sys, x_u, n_x, i_x_u, i_d_y, del_x_u, gradient_flag)

%JJ_LIN   Subsystem linearisation of a nonlinear ordinary differential equation system 
% 
%   JACO = JJ_LIN (SYS, X_U, N_X, I_X_U, I_D_Y)  
%   linearizes the system with the name SYS at the operating point
%   which is defined by the generalized input vector X_U.
%   N_X is the length of the original state vector X. 
%   It is needed for the disassembling of the X_U vector 
%   in the parameter list of the system calls.
%   The Jacobian matrix JACO only contains the subsystem defined by the
%   index vectors I_X_U and I_D_Y.
%
%   JACO = JJ_LIN (SYS, X_U, N_X, I_X_U, I_D_Y, DEL_X_U)
%   additionally allows the specification of the perturbation levels
%   DEL_X_U to be used for the gradient calculations.
%
%   JACO = JJ_LIN (SYS, X_U, N_X, I_X_U, I_D_Y, DEL_X_U, GRADIENT_FLAG)
%   additionally allows the specification of the difference stencil to use.
%   GRADIENT_FLAG can be defined as a vector of N = length(i_x_u) elements.
%   The values must be {-1,0,1}:
%   -1: backward differences,
%    0: central differences (default),
%    1: forward differences.
%   GRADIENT_FLAG can also be passed as scalar, in which case the same
%   difference stencil is applied for all trim variables.
%
%   Function dependencies: getGeneralizedOutputVector (for model evaluation)
%
%   Copyright 2000-2004, J. J. Buchholz, Hochschule Bremen, buchholz@hs-bremen.de


%   Version 1.3     22.07.2020, Daniel Kiehn, DLR-FT, daniel.kiehn@dlr.de
%   - added pre-allocation of the Jacobian matrix
%   - index vectors can now be column or row vectors
%   - added optional input gradient_flag
%   - Moved the model evaluation to new function getGeneralizedOutputVector


% If the user did not define the perturbation levels
if nargin < 6

	% Use default perturbation levels
	del_x_u = 1e-6*(1 + abs(x_u)); 
end

% If the user did not provide the gradient flag
if nargin < 7
	gradient_flag = zeros(size(i_x_u)); % Use central differences by default
else
	if isscalar(gradient_flag)
		gradient_flag = gradient_flag*ones(size(i_x_u));
	elseif length(gradient_flag) ~= length(i_x_u)
		errmsg = ['The length of gradient_flag (',num2str(length(gradient_flag)),')',...
				  ' does not match the length of i_x_u (',num2str(length(i_x_u)),').'];
		error(errmsg);
	end
end


% Determine vector lengths
n_i_x_u = length (i_x_u);
n_i_d_y = length (i_d_y);

% Set time to zero explicitly (assume a time-invariant system)
t = 0;

% Pre-allocate Jacobian matrix
jaco = zeros(n_i_d_y,n_i_x_u);

% Loop over all generalized inputs: generate one column of the Jacobian for
% every generalized input to be linearized and build up the matrix columnwise
for ii = 1:n_i_x_u

	% Save whole generalized input vector in x_u_left and x_u_right,
	% because we only want to waggle some specific generalized inputs
	x_u_left  = x_u;
	x_u_right = x_u;

	% Waggle the i-th generalized input
	x_u_left(i_x_u(ii))  = x_u(i_x_u(ii)) - del_x_u(i_x_u(ii));
	x_u_right(i_x_u(ii)) = x_u(i_x_u(ii)) + del_x_u(i_x_u(ii));

	% Central difference
	if gradient_flag(ii) == 0
		d_y_left  = getGeneralizedOutputVector(sys, x_u_left,  n_x, t);
		d_y_right = getGeneralizedOutputVector(sys, x_u_right, n_x, t);

		jaco_column = (d_y_right(i_d_y) - d_y_left(i_d_y))/(2*del_x_u(i_x_u(ii)));

	% Right-sided difference
	elseif gradient_flag(ii) == 1
		d_y_left  = getGeneralizedOutputVector(sys, x_u,       n_x, t);
		d_y_right = getGeneralizedOutputVector(sys, x_u_right, n_x, t);

		jaco_column = (d_y_right(i_d_y) - d_y_left(i_d_y))/(del_x_u(i_x_u(ii)));

	% Left-sided difference
	elseif gradient_flag(ii) == -1
		d_y_left  = getGeneralizedOutputVector(sys, x_u_left, n_x, t);
		d_y_right = getGeneralizedOutputVector(sys, x_u,      n_x, t);

		jaco_column = (d_y_right(i_d_y) - d_y_left(i_d_y))/(del_x_u(i_x_u(ii)));

	else
		errmsg = ['Wrong value of gradient_flag for generalized input #',...
				  num2str(ii),': ',num2str(gradient_flag(ii)),'. ',...
				  'Only {-1,0,1} are accepted!'];
		error(errmsg);
	end

	jaco(:,ii) = jaco_column;
end

end
