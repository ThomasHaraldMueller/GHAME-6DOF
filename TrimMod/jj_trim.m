function ...
  [x_tr, u_tr, d_tr, y_tr, varargout] = jj_trim (sys, ...
  x, u, d, y, ...
  i_x, i_u, i_d, i_y, ...
  x_nam, u_nam, d_nam, y_nam, ...
  del_x_max, del_u_max, del_x_lin, del_u_lin, ...
  options, ...
  x_min, x_max, u_min, u_max)

%JJ_TRIM   Trim point determination of a nonlinear ordinary differential
%equation system
%
%   [X_TR, U_TR, D_TR, Y_TR] = JJ_TRIM (SYS, X, U, D, Y, I_X, I_U, I_D, I_Y, X_NAM, U_NAM, D_NAM, Y_NAM)
%   trims the system SYS towards an operating point defined by the elements
%   of
%   X (state vector), U (input vector), D (derivative of the state vector),
%   and Y (output vector).
%   I_X and I_U define the indices of the so-called trim variables, which
%   are those states and inputs the trim algorithm has to find
%   the appropriate values for. Specified values of the trim variables
%   are taken as initial starting guesses for the iteration.
%   I_D and I_Y are the indices of the so-called trim requirements, which
%   the trim algorithm has to satisfy. The values of the other D(i) and Y(i)
%   do not matter.
%   X_NAM, U_NAM, D_NAM, and Y_NAM are cell arrays of the form
%   X_NAM = {'state_1'; 'state_2'; ...}
%   containing the names of the states, inputs, derivatives, and outputs.
%   The names can be chosen arbitrarily. They are used only to identify
%   linear dependent trim variables or trim requirements.
%
%   IMPORTANT: o There have to be as many trim variables as there are
%                trim requirements.
%
%              o All vectors (and cell arrays) have to be column vectors.
%
%   [X_TR, U_TR, D_TR, Y_TR, ERRFLG] = JJ_TRIM (SYS, ...)
%   ERRFLG = 0 if trimming was successful,
%   ERRFLG = 1 if an error occurred during trimming.
%              The returned trim result is invalid.
%
%   [X_TR, ..., ERRFLG, INFOSTRUCT] = JJ_TRIM (SYS, ...)
%   additionally returns an INFOSTRUCT structure containing
%   - x: the progression of all states,
%   - u: the progression of all inputs,
%   - d: the progression of all state derivatives,
%   - y: the progression of all outputs
%   over all iterations. They are of size N x (N_iter+1) where N is the
%   respective vector length (number of states, inputs etc.) and N_iter is
%   the number of iterations. The first column of each quantity represents
%   the initial point (before the first iteration).
%
%   - n_iter: the total number of iterations
%
%   - trimVars: the list of trim variables (cell array of strings)
%
%   - trimReqs: the list of trim requirements (cell array of strings)
%
%   - cost:   the progression of the cost function (length: N_iter+1)
%
%   - jaco:   the progression of the Jacobian as a cell array of N_iter
%			  cells, each cell containing the Jacobian (size N_t_v x N_t_v)
%             of the current iteration, where N_t_v is the number of trim
%			  variables / trim requirements.
%             The rows represent trim variables (as listed in trimVars),
%             the columns represent the trim requirements (as listed in
%             trimReqs).
%
%             Note that the Jacobian is not updated if the cost function
%             deteriorated from the previous step (StepType 2, see below).
%             In this case, an empty  matrix is stored instead.
%
%   - StepType: enumeration for the type of each step:
%         1: regular iteration
%         2: bisection to recover from a previous iteration step which lead
%            to an increased cost function
%         3: bisection to recover from a previous iteration step which lead
%            to a rank-deficient Jacobian even though the Jacobian had full
%            rank at the/a previous step
%    The i-th entry of StepType corresponds to the step from x(i) and u(i)
%    to x(i+1) and u(i+1).
%
%   - BisecCounter: number of bisections (including the current step) since
%                   the last successful iteration
%
%   - LimitedStep: flag indicating if the step size was truncated:
%     	 0: step was not truncated
%        1: step was truncated to the user-specified step size (DEL_X_MAX
%           and/or DEL_U_MAX)
%        2: step was truncated to avoid violation of upper or lower bounds
%           (X_MIN, X_MAX, U_MIN, or U_MAX)
%        3: step was truncated because of step size limitations AND hard
%           bounds.
%
%   - ExitCode: enumeration for the (reasons of) success or failure of the
%               trimming process:
%     	 0: successful trim
%        1: trim failed:  singular Jacobian at problem setup
%        2: trim failed:  max. number of bisections reached when trying to
%                         restore a regular Jacobian
%        3: trim failed:  max. number of bisections reached when trying to
%                         reduce the cost function
%        4: trim failed:  max. number of iterations reached
%        5: trim failed:  no more step possible in the required direction
%                         without violating the specified bounds
%
% 	- OutputMessage: string or cell of strings containing a detailed message
%                    about the trim result.
%
%
%   [X_TR, ...] = JJ_TRIM (..., Y_NAM, DEL_X_MAX, DEL_U_MAX)
%   allows the additional specification of maximum (absolute) step sizes of
%   states and inputs between consecutive iterations.
%   The lengths of DEL_X_MAX and DEL_U_MAX must be equal to those of
%   X and U respectively.
%   The default values of DEL_X_MAX and DEL_U_MAX are 1e42.
%
%   [X_TR, ...] = JJ_TRIM (..., DEL_U_MAX, DEL_X_LIN, DEL_U_LIN)
%   allows the additional specification of the state and input perturbations
%   to be used for the calculation of the Jacobian matrix (sensitivity matrix)
%   in the linearization procedure. Their length must match the respective
%   lengths of x and u. All elements must be nonzero and are automatically
%   treated as absolute values. The default values of DEL_X_LIN and DEL_U_LIN
%   are 1e-6*(1 + abs(x)) and 1e-6*(1 + abs(u)), respectively.
%
%   [X_TR, ...] = JJ_TRIM (..., DEL_U_LIN, OPTIONS)
%   allows the specification of additional options via designated fields of
%   the OPTIONS structure:
%   n_iter_max:   	the maximum number of iterations (default: 42)
%
%   cost_tbg:       the cost value to be gained (default: 1e-9)
%
%   CompileFlag:    flag {0,1,2} to control if/how to compile the model (default: 1)
%       o flag = 0: model is not compiled by jj_trim (hence needs to be precompiled)
%       o flag = 1: model is compiled by jj_trim using the 'compile' command.
%                   This is the default option in order to be backward
%                   compatible with jj_trim version 1.5 and lower.
%       o flag = 2: model is compiled by jj_trim using the 'lincompile'
%                   command. Recommended e.g. if the model incorporates time
%                   delay blocks to be used with direct feedthrough during
%                   linearization.
%
%   EnableMessages: flag to disable (0) or enable (~=0) message boxes and
%                   console messages. Messages are enabled by default. Note that
%                   disabling messages will make the function completely silent.
%
%	n_bisec_max:	Maximum number of consecutive bisections (default: 10)
%
%   If a field is not set (i.e. is not part of the structure), the
%   respective default value is used.
%   In order to maintain backward compatibility with jj_trim versions 1.5
%   and below, it is also possible to define OPTIONS as a 2-element vector.
%   In this case the following definition applies:
%	- OPTIONS(1): the maximum number of iterations (default: 42)
%   - OPTIONS(2): the cost value to be gained (default: 1e-9)
%   If OPTIONS is passed as a vector, it is not possible to set CompileFlag,
%   EnableMessages, or n_bisec_max, and their default values will be used.
%
%   OPTIONS can also be passed as an empty vector [], in which case the
%   default values will be used.
%
%
%   [X_TR, ...] = JJ_TRIM (..., OPTIONS, X_MIN, X_MAX, U_MIN, U_MAX)
%   allows the additional specification of fixed upper and lower boundaries
%   for states x and inputs u. They can be defined as vectors, in this case
%   their lengths must match the respective lengths of x and u.
%   You can use an empty vector [] to use the default values: -inf for
%   X_MIN and U_MIN, and inf for X_MAX and U_MAX, respectively.
%
%
%   Function dependencies:
%   1. jj_lin (for linearization)
%   2. getGeneralizedOutputVector (model evaluation)
%
%
%   Notes:
%   Unfortunately, Simulink claims the right to alter the order of the states
%   during the simulation (see https://mathworks.com/help/simulink/gui/initial-state.html).
%   MathWorks hence recommends not to use arrays for initialization. There
%   are different solutions to circumvent this; one solution is converting
%   the array to a structure (additionally features the corresponding state
%   names) and initializing the model using this structure as shown below:
%
%   u_tr_with_zero = [0, u_tr']; % add a zero as the time span trailing the input vector
%   x_tr_string = mat2str (x_tr, 42); % convert the state trim vector to a string with maximum precision:16
% 
%   set_param(sys,'LoadInitialState','on');  % enable the checkboxes in the Initial state...
%   set_param(sys,'LoadExternalInput','on'); % and Input menu entry of the current model
%   set_param(sys,'InitialState',x_tr_string ); % transfer the state trim vector to the corresponding edit text
%   set_param(sys,'SaveFormat','Structure'); % set the SaveFormat to Structure
% 
%   x_tr_structure = Simulink.BlockDiagram.getInitialState(sys); % reread the initial state back as a structure into a new variable
%   model_workspace = get_param(sys,'ModelWorkspace'); % retrieve a handle to the model workspace
% 
%   assignin(model_workspace,'initial_state_in_model_workspace',x_tr_structure); % write the initial state (as a structure)...
%   assignin(model_workspace,'initial_input_in_model_workspace',u_tr_with_zero); % and the initial input (as a vector) to the model
% 
%   set_param(app.model_name,'ExternalInput','initial_input_in_model_workspace'); % use these model workspace variables...
%   set_param(app.model_name,'InitialState','initial_state_in_model_workspace');  % in the Input and Initial state edit texts



%   Authors:
%   J. J. Buchholz, Hochschule Bremen,              buchholz@hs-bremen.de
%   W. Moennich,  	German Aerospace Center (DLR),	wulf.moennich@dlr.de
%   D. Kiehn,       German Aerospace Center (DLR),	daniel.kiehn@dlr.de



%   Version 2.3.1   16.01.2024, D. Kiehn
%   - After setting 'InitInArrayFormatMsg' to 'none' (in newer Matlab/Simulink
%     versions), the model is now marked as unchanged ('dirty' = 'off')
%     which causes Simulink not to ask whether or not the model is to be
%     saved, as was previously the case.


%   Version 2.3     03.01.2024, D. Kiehn
%   - Added the names of trim variables and trim requirements to the
%     infostruct output.


%   Version 2.2     25.02.2022, D. Kiehn
%   - Fixed a bug which caused the model not to be terminated in case the
%     trim failed due to hard bounds (i.e., exit code 5)
%   - It is now verified that the target cost function (cost_tbg) is non-negative.
%   - Added a note about the risks of directly using the trim result (as an array)
%     for the initialization of a model


%   Version 2.1     05.08.2021, D. Kiehn
%   - Fixed a bug which caused the step direction to change sign in case
%     the step had to be truncated to avoid violating the lower bounds
%     (X_MIN or U_MIN).
%   - OPTIONS can now be passed as an empty element.


%   Version 2.0     22.07.2020, D. Kiehn
%   - Individual checking of max. iteration step sizes and max. perturbations
%     for linearization.
%   - The suppression of the array initialization warning is now done via
%     try/catch to ensure compatibility with older Matlab versions.
%   - Added automatic bisection if the sensitivity matrix loses full rank
%     (but was initialized with full rank) - useful in case a saturation is
%     reached during the trimming process.
%   - The OPTIONS input is now a structure which can be used to
%       - control if/how to compile the Simulink model, granting the
%         ability to use precompiled models,
%       - define the maximum number of iterations and bisections,
%       - enable/disable all messages and warnings.
%   - Added optional inputs to specify fixed upper and lower boundaries for
%     states and inputs: X_MIN, X_MAX, U_MIN, U_MAX.
%   - Added optional output INFOSTRUCT with details about the trimming
%     process (progression of state/input vectors, cost function, Jacobian,
%     etc. over all iterations).
%   - Introduced a "fail-fast" philosophy: the function now returns an
%     error if it is called with nonsensical inputs, e.g. if the number of
%     trim variables doesn't match the number of trim requirements.
%   - Moved the model evaluation to new function getGeneralizedOutputVector
%
%   Note: default options of version 2.0 are those of jj_trim version 1.5
%         to maintain backward compatibility.


%   Version 1.5     07.11.2011, J. J. Buchholz
%
%   Array initialization warning suppressed.
%   Minor code optimizations.


%   Version 1.2.4   15.10.2008, W. Moennich
%
%   Comment added:
%   % The functionalities used here ('compile', 'term', ...)
%   % are documented in the Simulink manual under Simulation Commands/model
%   % and in Help(doc) unter Simulink/Functions/Simulation/model.


%   Version 1.2.3   20.05.2003, W. Moennich
%
%   Optional error flag as 5th output.


%   Version 1.2.2   26.09.2001, W. Moennich
%
%   Empty vectors now possible for DEL_X_MAX, DEL_U_MAX, DEL_X_LIN,
%   DEL_U_LIN.


%   Version 1.2.1   29.08.2001, W. Moennich
%
%   Check for maximum number of iterations moved inside the loop.


%   Version 1.2 	26.05.2000, J. J. Buchholz
%
%   The names of inputs, ... outputs for the error messages
%   are now transferred via the parameter list.
%
%   The precompilation and the release of the system is now done in JJ_TRIM.


%   Version 1.0     1998,       J. J. Buchholz
%
%   Initial version.



%% Initialize

% Open the system block diagram internally without bringing the system
% block diagram to front
dummy = eval (sys);

% Suppress the array initialization warning (only required in newer Matlab
% versions and will crash in older Matlab versions, hence the try/catch)
try
	set_param(sys, 'InitInArrayFormatMsg', 'None');
    set_param(sys, 'Dirty', 'off') % Otherwise, Simulink will treat the model as "changed"
catch
end

% Feed through all initial vectors, useful in case of an emergency exit
x_tr = x;           % States
u_tr = u;           % Inputs
d_tr = d;           % State derivatives
y_tr = y;           % Outputs
varargout{1} = [];  % Error flag
varargout{2} = {};  % Infostruct

% Determine lengths of basic vectors
n_x  = length (x);
n_u  = length (u);
n_d  = length (d);
n_y  = length (y);

n_i_x = length (i_x); % Number of trim variables (states)
n_i_u = length (i_u); % Number of trim variables (inputs)
n_i_d = length (i_d); % Number of trim requirements (state derivatives)
n_i_y = length (i_y); % Number of trim requirements (outputs)
%disp(x)
%disp(u)
% Assemble generalized input vector and generalized output vector
x_u = [x; u];
d_y = [d; y];
%disp(d_y)

% Determine length of generalized input vector
n_x_u = n_x + n_u;

% Assemble trim variable index vector and trim requirement index vector
i_t_v = [i_x; i_u + n_x];
i_t_r = [i_d; i_y + n_d];

% Determine number of trim variables and trim requirements
n_t_v = n_i_x + n_i_u;
n_t_r = n_i_y + n_i_d;

% Set evaluation time to zero explicitly (assume a time invariant system)
t = 0;


%% Read and check inputs

% Options have to be checked first to decide whether or not to display
% error messages
% If no options have been defined,
if nargin < 18

	% Set defaults which match jj_trim 1.5 and lower
	n_iter			= 52;	% Max. number of iterations (including bisection steps)
	cost_tbg		= 1e-0;	% Cost function to be gained
	CompileFlag		= 1;	% Controls if/how to compile the Simulink model (see function help)
	EnableMessages	= 1;	% Disables (0) or enables (~=0) pop-up message boxes
	N_bisec_max     = 10;	% Max. number of consecutive bisections
else
	if isempty(options) % if options are passed as empty vector, use defaults
		n_iter = 42;
		cost_tbg = 1e-9;
		CompileFlag = 1;
		EnableMessages = 1;
		N_bisec_max = 10;

	elseif isstruct(options)
		% If the options are defined as a struct, read the specified options
		% individually

		% Read max. number of iterations
		if isfield(options,'n_iter_max')
			n_iter = options.n_iter_max;
		else
			n_iter = 452;
		end

		% Read cost value to be gained
		if isfield(options,'cost_tbg')
			if options.cost_tbg >= 0
				cost_tbg = options.cost_tbg;
			else
				errordlg ('The cost value to be gained must not be negative!','Error');
				return
			end
		else
			cost_tbg = 1e-10;
		end

		% Read the compile flag
		if isfield(options,'CompileFlag')
			CompileFlag = options.CompileFlag;
		else
			CompileFlag = 1;
		end

		% Read the flag that enables/disables messages
		if isfield(options,'EnableMessages')
			EnableMessages = options.EnableMessages;
		else
			EnableMessages = 1;
		end

		% Read the max. number of consecutive bisections
		if isfield(options,'n_bisec_max')
			N_bisec_max = options.n_bisec_max;
		else
			N_bisec_max = 10;
		end

	else
		% If the "options" input is not a struct, try to access the vector
		% entries. This is required to maintain backward compatibility
		% with versions 1.5 and below.
		n_iter   = options(1);
		cost_tbg = options(2);
		EnableMessages = 1;
		CompileFlag	   = 1;
		N_bisec_max = 10;
	end
end


% Check if there are as many trim variables as there are trim requirements
if n_t_r ~= n_t_v

	l1 = ['The number of trim variables: ', int2str(n_t_v)];
	l2 = 'does not equal';
	l3 = ['the number of trim requirements: ', int2str(n_t_r)];
	l4 = ' ';
	l5 = 'The returned trim point is not valid.';

	outputMessage = {l1; l2; l3; l4; l5};
	errordlg (outputMessage,'Error');
	return
end


% If there are no trim variables and requirements
if ~n_t_r
	if EnableMessages
		% Prepare output message
		outputMessage = {'There should be at least';
						 'one trim variable and';
						 'one trim requirement.'};
		warndlg(outputMessage,'Nothing to trim');
	end
	return
end


% If no maximum step sizes have been defined for x and u
if nargin < 14

	% set defaults
	del_max = 1.e42*ones (n_x_u, 1);

	% If the maximum step sizes have been defined at least for x
else

	% Check del_x_max, set to default if empty
	if isempty(del_x_max)
		del_x_max = 1.e42*ones (n_x,1);
	elseif numel(del_x_max) ~= n_x
		outputMessage = {...
			['The number of elements in del_x_max (',num2str(numel(del_x_max)),')'];
			['must match the length of x (' num2str(n_x),').']};
		errordlg(outputMessage,'Error');
		return
	end

	% If no maximum step sizes have been defined for u,
	if nargin < 15
		del_u_max = 1.e42*ones (n_u, 1);
	else
		% Check del_u_max, set to default if empty
		if isempty(del_u_max)
			del_u_max = 1.e42*ones (n_u, 1);
		elseif numel(del_u_max) ~= n_u
			outputMessage = {...
				['The number of elements in del_u_max (',num2str(numel(del_u_max)),')'];
				['must match the length of u (' num2str(n_u),').']};
			errordlg (outputMessage,'Error');
			return
		end
	end

	% Assemble generalized maximum step vector
	del_max = [del_x_max(:); del_u_max(:)];
end


% If no step sizes for the linearization have been defined,
if nargin < 16
	del_lin = 1e-6*(1 + abs(x_u)); % set defaults
else
	% Check del_x_lin, set to default if empty
	if isempty(del_x_lin)
		del_x_lin = 1e-6*(1 + abs(x));
	elseif numel(del_x_lin) ~= n_x
		outputMessage = {...
			['The number of elements in del_x_lin (',num2str(numel(del_x_lin)),')'];
			['does not match the length of x (' num2str(n_x),').']};
		errordlg (outputMessage,'Error');
		return
	end

	% If no maximum step sizes have been defined for u,
	if nargin < 17
		del_u_lin = 1e-6*(1 + abs(u)); % set defaults
	else
		% Check del_u_lin, set to default if empty
		if isempty(del_u_lin)
			del_u_lin = 1e-6*(1 + abs(u));
		elseif numel(del_u_lin) ~= n_u
			outputMessage = {...
				['The number of elements in del_u_lin (',num2str(numel(del_u_lin)),')'];
				['does not match the length of u (' num2str(n_u),').']};
			errordlg(outputMessage,'Error');
			return
		end
	end

	% Assemble generalized perturbation vector for linearization
	% Make sure only positive values are used
	del_lin = [abs(del_x_lin(:)); abs(del_u_lin(:))];
end


% Check upper and lower bounds for x and u
if nargin < 19

	% If no bounds were specified at all, set defaults
	x_min = -inf*ones(n_x,1);
	x_max =  inf*ones(n_x,1);
	u_min = -inf*ones(n_u,1);
	u_max =  inf*ones(n_u,1);
else
	% If at least x_min was specified,
	if isempty(x_min)
		x_min = -inf*ones(n_x,1);
	elseif numel(x_min) == n_x
		x_min = x_min(:);
	else
		outputMessage = {...
			['The number of elements in x_min (',num2str(numel(x_min)),')'];
			['does not match the length of x (' num2str(n_x),').']};
		errordlg(outputMessage,'Error');
		return
	end

	% Check conformance of x_max with initial x
	if any(x(:) < x_min(:))
		errordlg('The initial state vector x violates the bounds of x_min!','Error');
		return
	end
end


% Check x_max
if nargin > 19
	if isempty(x_max)
		x_max = inf*ones(n_x,1);
	elseif numel(x_max) == n_x
		x_max = x_max(:);
	else
		outputMessage = {...
			['The number of elements in x_max (', num2str(numel(x_max)),')'];
			['does not match the length of x (' num2str(n_x),').']};
		errordlg(outputMessage,'Error');
		return
	end

	% Check conformance of x_max with initial x
	if any(x(:) > x_max(:))
		errordlg('The initial state vector x violates the bounds of x_max!','Error');
		return
	end
else
	x_max =  inf*ones(n_x,1);
end


% Check u_min
if nargin > 20
	if isempty(u_min)
		u_min = -inf*ones(n_u,1);
	elseif numel(u_min) == n_u
		u_min = u_min(:);
	else
		outputMessage = {...
			['The number of elements in u_min (',num2str(numel(u_min)),')'];
			['does not match the length of u (' num2str(n_u),').']};
		errordlg(outputMessage,'Error');
		return
	end

	% Check conformance of u_min with initial u
	if any(u < u_min)
		errordlg('The initial input vector u violates the bounds of u_min!','Error');
		return
	end
else
	u_min = -inf*ones(n_u,1);
end


% Check u_max
if nargin > 21
	if isempty(u_max)
		u_max = inf*ones(n_u,1);
	elseif numel(u_max) == n_u
		u_max = u_max(:);
	else
		outputMessage = {...
			['The number of elements in u_max (', num2str(numel(u_max)),')'];
			['does not match the length of u (' num2str(n_u),').']};
		errordlg(outputMessage,'Error');
		return
	end

	% Check conformance of u_max with initial u
	if any(u > u_max)
		errordlg('The initial input vector u violates the bounds of u_max!','Error');
		return
	end
else
	u_max = inf*ones(n_u,1);
end


x_u_min = [x_min; u_min];
x_u_max = [x_max; u_max];


%% Core algorithm

% Pre-allocate infostruct (only required if 6th output is desired)
if nargout > 5
	infostruct.x    = zeros(n_x,n_iter+1);
	infostruct.u    = zeros(n_u,n_iter+1);
	infostruct.d    = zeros(n_d,n_iter+1);
	infostruct.y    = zeros(n_y,n_iter+1);
	infostruct.cost = zeros(1,n_iter+1);

	infostruct.jaco = cell(1,n_iter);

	infostruct.StepType     = zeros(1,n_iter);
	infostruct.BisecCounter = zeros(1,n_iter);
	infostruct.LimitedStep  = zeros(1,n_iter);

	infostruct.n_iter   = 0;
	infostruct.ExitCode = 0;
    
    % Additionally, store the names of trim variables and requirements
    x_u_nam = [x_nam; u_nam];               % Names of all generalized inputs
    d_y_nam = [d_nam; y_nam];               % Names of all generalized outputs
    infostruct.trimVars = x_u_nam(i_t_v);	% Names of trim variables
    infostruct.trimReqs = d_y_nam(i_t_r); 	% Names of trim requirements
end

% Set old cost value to infinity, in oder to definitely have 
% an improvement with the first try
cost_old = inf;

% Precompile system.
% Unfortunately, this is necessary because only precompiled systems can be evaluated.
% If the trim algorithm is aborted without the corresponding "Release system" command
% the next precompilation attempt will lead to an error and the simulation cannot
% be started.
% The system then has to be released manually (maybe more than once!) with:
% model_name ([], [], [], 'term')
% The functionalities used here ('compile', 'term', ...) 
% are documented in the Simulink manual under Simulation Commands/model
% and in Help(doc) unter Simulink/Functions/Simulation/model.
if CompileFlag == 1
	feval (sys, [], [], [], 'compile');
elseif CompileFlag == 2
	feval (sys, [], [], [], 'lincompile');
end

% Initialize bisection counter
i_bisec = 0;

% Loop over maximum n_iter iterations
for i_iter = 0 : n_iter

	% Calculate generalized output vector at the current trim point.
	d_y_tr = getGeneralizedOutputVector(sys,x_u,n_x,t);
    %disp(d_y_tr)

	% Update function outputs
	x_tr = x_u(1:n_x,1);
	u_tr = x_u(n_x+1:end,1);
	d_tr = d_y_tr(1:n_d,1);
	y_tr = d_y_tr(n_d+1:end,1);

	% Calculate differences between required and current generalized output vectors
	del_d_y = d_y - d_y_tr;

	% Pick trim requirements only
	del_t_r = del_d_y(i_t_r);

	% Cost value is the maximum element of the trim requirement error vector
	cost = max(abs(del_t_r));
    
    % Check for NaNs in the cost function
    if any(isnan(del_t_r)) && EnableMessages
        disp(['Warning: NaN detected in trim requirements in iteration ',num2str(i_iter)])
    end

	% Update infostruct (if output is desired)
	if nargout > 5
		infostruct.x(1:n_x,i_iter+1) = x_u(1:n_x,1);
		infostruct.u(1:n_u,i_iter+1) = x_u(n_x+1:end,1);
		infostruct.d(1:n_d,i_iter+1) = d_y_tr(1:n_d,1);
		infostruct.y(1:n_y,i_iter+1) = d_y_tr(n_d+1:end,1);

		infostruct.cost(1,i_iter+1) = cost;
	end

	% If current cost value has become smaller 
	% than the cost value to be gained
	if cost < cost_tbg

		% Assemble output message
		l1 = ['A cost value of ',num2str(cost)];
		l2 = ['has been gained after ', int2str(i_iter), ' iteration(s).'];
		h1 = 'Success';

		outputMessage = [l1,' ',l2];

		% Display success message in console
		if EnableMessages
            helpdlg({l1, l2}, h1);
		end

		% Release system
		if (CompileFlag == 1) || (CompileFlag == 2)
			feval (sys, [], [], [], 'term');
		end

		% If a 5th output parameter is used,
		if nargout > 4

			% Return errorflag as 5th parameter
			errflg = 0;
			varargout{1} = errflg;

			% If a 6th output parameter is used,
			if nargout > 5

				% Update infostruct and remove unassigned parts
				infostruct.ExitCode	= 0; % Exit code 0: successful trim
				infostruct.OutputMessage = outputMessage;
				infostruct = truncateInfoStruct(infostruct,i_iter);

				% Assign infostruct to second output
				varargout{2} = infostruct;
			end
		end

		% Trimming was successful
		return
	end

	% If the maximum number of iterations is reached,
	% output error message and abort program
	if i_iter == n_iter

		l1 = ['Maximum number of iterations reached: ', int2str(i_iter)];
		l2 = 'Program was aborted';
		l3 = ['with a cost value of: ', num2str(cost)];

		outputMessage = {l1; l2; l3};

		if EnableMessages
			errordlg (outputMessage, 'Program aborted');
		end

		% Release system
		if (CompileFlag == 1) || (CompileFlag == 2)
			feval (sys, [], [], [], 'term');
		end

		% If a 5th output parameter is used,
		if nargout > 4

			% Update error flag and return it as 5th parameter
			errflg = 1;
			varargout{1} = errflg;

			% If a 6th output parameter is used,
			if nargout > 5

				% Update infostruct and remove unnecessary fields
				infostruct.ExitCode	= 4; % Exit code 4: max. iterations reached
				infostruct.OutputMessage = outputMessage;
				infostruct = truncateInfoStruct(infostruct,i_iter);

				% Assign infostruct to second (optional) output
				varargout{2} = infostruct;
			end
		end

		% Game over: max. iterations reached
		return
	end


	% If an improvement has been obtained
	% with respect to the last point,
	if cost < cost_old

		% Assemble the gradient flag vector. This vector stores information
		% about the difference stencils to use during linearization:
		% -1 for backward difference,
		%  0 for central difference (default),
		% +1 for forward difference.
		% This way it is ensured that the evaluations during linearization
		% do not violate the bounds of x_u_min and x_u_max.
		gradient_flag = assembleGradientFlags(x_u,i_t_v,del_lin,x_u_min,x_u_max);

		% Linearize relevant subsystem at current operating point
		jaco = jj_lin (sys, x_u, n_x, i_t_v, i_t_r, del_lin, gradient_flag);

		% Assign current Jacobian to infostruct
		if nargout > 5
			infostruct.jaco{i_iter+1} = jaco;
		end

		% Singular Value Decomposition of the sensitivity matrix
		[u, s, v] = svd (jaco);

		% A singular value is assumed to be "zero", if it is 1e12 times smaller 
		% than the maximum singular value. Such a singular value indicates a rank deficiency.
		sv_min = s(1,1)*1e-12;

		% Find the indices of those "zero-singular-values"
		i_sv_zero = find (abs (diag (s)) <= sv_min);


		% If there are no zero-singular-values,
		if isempty(i_sv_zero)

			% accept and save this new point.
			% Important for a possible step size bisection later on
			x_u_old = x_u;

			% Save the cost value of this new point for a comparison later on
			cost_old = cost;

			% Reset step size bisection counter
			i_bisec = 0;

			% Assuming a linear system, the variation of the trim variables
			% necessary to compensate the trim requirements error can directly
			% be calculated by the inversion of the linear subsystem model
			% (differential equations and output equations)
			del_t_v = jaco\del_t_r;

			[del_t_v,stepLengthViolated] = forceConformanceWithMaximumStepLength(del_t_v,del_max(i_t_v));

			[del_t_v,hardBoundsViolated] = forceConformanceWithHardLimits(x_u,del_t_v,i_t_v,x_u_min,x_u_max);            

			% Set the LimitedStep flag accordingly
			if stepLengthViolated && hardBoundsViolated
				LimitedStep = 3;
			elseif hardBoundsViolated
				LimitedStep = 2;
			elseif stepLengthViolated
				LimitedStep = 1;
			else
				LimitedStep = 0;
			end

			% If at least one element of the step is nonzero
			if ~all(del_t_v == 0)

				% Accept the step and calculate new trim point.
				% Always use old value *before* first bisection
				x_u(i_t_v) = x_u_old(i_t_v) + del_t_v;

				% Disassemble the generalized input vector
				x_tr = x_u(1:n_x);
				u_tr = x_u(n_x+1:end);
				d_tr = d_y_tr(1:n_d);
				y_tr = d_y_tr(n_d+1:end);

				% Update step type: 1 = regular iteration
				if nargout > 5
					steptype = 1;
					infostruct.StepType(i_iter+1)       = steptype;
					infostruct.BisecCounter(1,i_iter+1) = i_bisec;
					infostruct.LimitedStep(1,i_iter+1)  = LimitedStep;
				end
			else
				% Game over: no step in the required direction possible (dead end)

				% Generate output message
				l1 = ['Trim failed in iteration ', num2str(i_iter),':'];
				l2 = 'No step in the required direction allowed because';
				l3 = 'the calculated step would violate the specified bounds.';
				l4 = ' ';
				l5 = 'Either the trim requirements cannot be met';
				l6 = 'or the specified min./max. bounds are too strict.';
				l7 = ' ';
				l8 = 'Try trimming without hard bounds';
				l9 = 'or try different initial values.';

				outputMessage = {l1; l2; l3; l4; l5; l6; l7; l8; l9};

				% Display error message if desired
				if EnableMessages
					errordlg (outputMessage,'Dead end');
				end

				% Release system
				if (CompileFlag == 1) || (CompileFlag == 2)
					feval (sys, [], [], [], 'term');
				end

				if nargout > 4
					errflg = 1;
					varargout{1} = errflg;

					if nargout > 5

						infostruct.ExitCode = 5; % Exit code 5: dead-end
						infostruct.OutputMessage = outputMessage;
						infostruct = truncateInfoStruct(infostruct,i_iter);

						varargout{2} = infostruct;
						return
					end
				end
			end
		else % if the Jacobian is singular

			% Update step type and bisection counter
			if nargout > 5
				steptype = 3;
				infostruct.StepType(i_iter+1)       = steptype;
				infostruct.BisecCounter(1,i_iter+1) = i_bisec+1;
				infostruct.LimitedStep(1,i_iter+1)  = 0;
			end


			% If the initial Jacobian is the problem
			if i_iter == 0

				% Generate/display output message
				outputMessage = generateErrorMessage(x_nam,u_nam,d_nam,y_nam,u,v,i_t_v,i_t_r,sv_min,i_sv_zero,0,EnableMessages);

				% Update optional outputs and return
				if nargout > 4

					% Update error flag and return it as 5th parameter
					errflg = 1;
					varargout{1} = errflg;

					% If a 6th output parameter is used,
					if nargout > 5

						% Update infostruct and remove unassigned parts
						infostruct.ExitCode	= 1; % Exit code 1: singular initial Jacobian
						infostruct.OutputMessage = outputMessage;
						infostruct = truncateInfoStruct(infostruct,i_iter);

						% Store Jacobian (necessary because i_iter = 0)
						infostruct.jaco{1} = jaco;

						% Assign infostruct to second (optional) output
						varargout{2} = infostruct;
					end
				end

				% Release system
				if (CompileFlag == 1) || (CompileFlag == 2)
					feval (sys, [], [], [], 'term');
				end

				% Game over: singular Jacobian at initial point
				return
			end


			% If step size has not been bisected N_bisec_max times yet,
			if i_bisec < N_bisec_max

				% bisect step size and change sign
				del_t_v = del_t_v/2;

				% and increment bisection counter
				i_bisec = i_bisec + 1;

				% Calculate new trim point.
				% Always use last value where Jacobian was regular and cost
				% had decreased
				x_u(i_t_v) = x_u_old(i_t_v) + del_t_v;

			% If step size has already been bisected N_bisec_max times before,
			else
				% Prepare error message and abort

				% Assemble output message
				l1 = ['Maximum number of consecutive bisections (',num2str(N_bisec_max),') reached'];
				l2 = 'while trying to restore a regular Jacobian.';
				l3 = ['Program was aborted after ', int2str(i_iter), ' iteration(s)'];
				l4 = ['with a cost value of ', num2str(cost), '.'];
				l5 = 'Try different initial values.';
				l6 = 'Or try to reduce step sizes.';

				outputMessage = {l1; l2; l3; l4; l5; l6};

				% Display error message if desired
				if EnableMessages
					errordlg (outputMessage,'Program aborted');
				end

				% Release system
				if (CompileFlag == 1) || (CompileFlag == 2)
					feval (sys, [], [], [], 'term');
				end

				% Game over: max. number of bisections reached
				if nargout > 4
					errflg = 1;
					varargout{1} = errflg;

					if nargout > 5

						% Exit code 2: max. number of bisections reached
						% when trying to restore a regular Jacobian
						% Update infostruct and remove unassigned parts
						infostruct.ExitCode	= 2;
						infostruct.OutputMessage = outputMessage;
						infostruct = truncateInfoStruct(infostruct,i_iter);

						varargout{2} = infostruct;
					end
				end

				return
			end
		end

	% If the cost function has degraded
	else

		% Store step type 2 and bisection counter
		if nargout > 5
			steptype = 2;
			infostruct.StepType(i_iter+1)       = steptype;
			infostruct.BisecCounter(1,i_iter+1) = i_bisec+1;
			infostruct.LimitedStep(1,i_iter+1)  = 0;

			% Store Jacobian as empty matrix to make clear that it has not
			% been updated since last successful iteration step
			infostruct.jaco{i_iter+1} = [];
		end

		if i_bisec < N_bisec_max

			% bisect step size and change sign
			del_t_v = del_t_v/2;

			% and increment bisection counter
			i_bisec = i_bisec + 1;

			% Calculate new trim point.
			% Always use last value where Jacobian was regular and cost
			% had decreased
			x_u(i_t_v) = x_u_old(i_t_v) + del_t_v;
		else

			% Prepare error message and abort

			% Assemble output message
			l1 = ['Maximum number of consecutive bisections (',num2str(N_bisec_max),') reached'];
			l2 = 'while trying to reduce the cost function.';
			l3 = ['Program was aborted after ', int2str(i_iter), ' iteration(s)'];
			l4 = ['with a cost value of ', num2str(cost), '.'];
			l5 = 'Try different initial values.';
			l6 = 'Or try to reduce step sizes.';

			outputMessage = {l1; l2; l3; l4; l5; l6};

			% Display error message if desired
			if EnableMessages
				errordlg (outputMessage,'Program aborted');
			end

			% Release system
			if (CompileFlag == 1) || (CompileFlag == 2)
				feval (sys, [], [], [], 'term');
			end

			% Game over: max. number of bisections reached
			if nargout > 4
				errflg = 1;
				varargout{1} = errflg;

				if nargout > 5

					% Exit code 3: max. number of bisections reached while
					% trying to reduce the cost function
					infostruct.ExitCode	= 3;
					infostruct.OutputMessage = outputMessage;
					infostruct = truncateInfoStruct(infostruct,i_iter);

					varargout{2} = infostruct;
				end
			end

			return
		end
	end
end

end



function [del_t_v,stepLengthViolated] = forceConformanceWithMaximumStepLength(del_t_v,del_t_v_max)
% This function ensures that the calculated step length del_t_v does not
% violate the maximum step length del_t_v_max.

	% Calculate maximum ratio between necessary and allowed trim
	% step sizes
	ratio_t_v = del_t_v ./ del_t_v_max;
	max_rat = max(abs(ratio_t_v));

	% If allowed step size has been exceeded,
	if max_rat > 1
		% scale all state and input step sizes, 
		% in order to exploit most of the allowed step size
		del_t_v = del_t_v/max_rat;

		stepLengthViolated = 1;
	else
		stepLengthViolated = 0;
	end
end



function [del_t_v,hardBoundsViolated] = forceConformanceWithHardLimits(x_u,del_t_v,i_t_v,x_u_min,x_u_max)
% This function ensures that the calculated step length del_t_v for the
% current x_u does not violate the bounds of x_u_min and x_u_max.

	hardBoundsViolated = 0;

	ratio_t_v = del_t_v./(x_u_max(i_t_v) - x_u(i_t_v));
	max_rat = max(ratio_t_v);

	if max_rat > 1
		del_t_v = del_t_v/max_rat;
		hardBoundsViolated = 1;
	end

    
	ratio_t_v = del_t_v./(x_u(i_t_v) - x_u_min(i_t_v));
	min_rat = min(ratio_t_v);
 
	if min_rat < -1
		del_t_v = del_t_v/-min_rat;
		hardBoundsViolated = 1;
	end
end



function gradient_flag = assembleGradientFlags(x_u,i_t_v,del_lin,x_u_min,x_u_max)
% This function assembles the vector which declares if the gradients of a
% trim variable should be calculated using central difference (flag = 0),
% forward difference (flag = 1) or backward difference (flag = -1). Central
% differences are used by default, and forward/backward differences are
% only used if the linearization step from the current working point would
% violate the upper or lower bounds x_u_max/x_u_min.
% Note: The implemented logic only works if all elements of the linearization
%       step vector del_lin are positive.

	n_tv = length(i_t_v);

	gradient_flag = zeros(n_tv,1);

	for ii = 1:n_tv

		if (x_u(i_t_v(ii)) + del_lin(i_t_v(ii))) > x_u_max(i_t_v(ii))
			gradient_flag(ii) = -1;
		elseif (x_u(i_t_v(ii)) - del_lin(i_t_v(ii))) < x_u_min(i_t_v(ii))
			gradient_flag(ii) = +1;
		else
			gradient_flag(ii) = 0;
		end
	end
end



function infostruct = truncateInfoStruct(infostruct_in,n_iter)
% This function truncates the fields of infostruct to the number of
% iterations actually performed. This is necessary since it was
% pre-allocated with a maximum size, but the final structure should
% contain only the actual iterations (n_iter plus the initial point).

	% Copy/paste structure - some entries will remain the same
	infostruct = infostruct_in;

	infostruct.n_iter = n_iter;

	infostruct.x = infostruct_in.x(:,1:n_iter+1);
	infostruct.u = infostruct_in.u(:,1:n_iter+1);
	infostruct.d = infostruct_in.d(:,1:n_iter+1);
	infostruct.y = infostruct_in.y(:,1:n_iter+1);

	infostruct.cost = infostruct_in.cost(1:n_iter+1);

	% Truncate step information
	infostruct.StepType		= infostruct_in.StepType(1:n_iter);
	infostruct.BisecCounter = infostruct_in.BisecCounter(1:n_iter);
	infostruct.LimitedStep  = infostruct_in.LimitedStep(1:n_iter);

	infostruct.jaco = infostruct_in.jaco(1:n_iter);
end



function errorMessage = generateErrorMessage(x_nam,u_nam,d_nam,y_nam,u,v,i_t_v,i_t_r,sv_min,i_sv_zero,bisectionWasUsed,enableMessages)
% This function generates the message that is displayed in case of a singular
% Jacobian matrix. It informs the user which trim variables and/or trim
% requirements might be the problem.

	% Generate top and bottom part of the error message to distinguish
	% whether or not bisection was used to restore a regular Jacobian
	% Note: the middle part of the message is common in both cases.
	if bisectionWasUsed == 0

		messageTitle = 'Singular initial Jacobian matrix!';

		upperPart = ' ';

		lowerPart = {...
			' ';
			'Choose different trim variables and/or trim requirements.';
			'Or try different initial values.';
			' ';
			'The returned trim point is not valid.';
			'If you are using trimmod, you can use the';
			'Untrim menu entry to return to the pre-trim state.'};

	else

		messageTitle = 'Singular Jacobian matrix!';

		upperPart = {...
			'Originally the trim problem was correctly';
			'set up (regular Jacobian matrix).';
			'This property was lost along the way and';
			'could not be restored by bisection.'};

		lowerPart = {...
			' ';
			'The returned trim point is not valid.';
			'If you are using trimmod, you can use the';
			'Untrim menu entry to return to the pre-trim state.'};
	end

	% Initialize the core part of the error message
	middlePart = [];

	% Assemble cell arrays containing the names of all trim variables and trim requirements
	trim_variables = [x_nam; u_nam];
	trim_requirements = [d_nam; y_nam];

	% Loop over all zero-singular-values
	for i_sv = i_sv_zero'

		% Find those elements of the corresponding singular vectors that are not "zero"
		u_sing = find (abs (u(:,i_sv)) > sv_min);
		v_sing = find (abs (v(:,i_sv)) > sv_min);

		% Separating empty line
		l10 = {' '};

		% If there is only one zero element in the left singular vector,
		if length (u_sing) == 1

			% prepare the corresponding error message
			l11 = {'The trim requirement'};
			l12 = trim_requirements(i_t_r(u_sing));
			l13 = {'could not be affected by any trim variable.'};

		% If there are more than one zero element in the left singular vector
		else

			% prepare the corresponding error message
			l11 = {'The trim requirements'};
			l12 = trim_requirements(i_t_r(u_sing));
			l13 = {'linearly depend on each other.'};
		end

		% Separating empty line
		l14 = {' '};

		% If there is only one zero element in the right singular vector,
		if length (v_sing) == 1

			% prepare the corresponding error message
			l15 = {'The trim variable'};
			l16 = trim_variables(i_t_v(v_sing));
			l17 = {'does not affect any trim requirement.'};

		% If there are more than one zero element in the right singular vector
		else

			% prepare the corresponding error message
			l15 = {'The trim variables'};
			l16 = trim_variables(i_t_v(v_sing));
			l17 = {'linearly depend on each other.'};
		end

		% Display error message
		if enableMessages
			errordlg([upperPart; l10; l11; l12; l13; l14; l15; l16; l17; lowerPart],messageTitle);
		end

		% Enhance the middle part of the complete error message
		middlePart = [middlePart; l10; l11; l12; l13; l14; l15; l16; l17;]; %#ok<AGROW>
	end

	% Assemble final output message
	errorMessage = [messageTitle; upperPart; middlePart; lowerPart];
end

