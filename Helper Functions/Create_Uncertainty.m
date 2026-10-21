function UNC = Create_Uncertainty(Data, scope, level, seed)
% CREATE_UNCERTAINTY  Perturbed aerodynamic database for one Monte Carlo run,
%   with a selectable scope of uncertain coefficients.
%
%   UNC = Create_Uncertainty(Data, scope, level, seed)
%
%   Each selected aerodynamic table is scaled by one uniformly distributed
%   factor (1 + delta), delta in [-level, level]. One factor per table keeps
%   the shape of the table, so slopes and the derivatives stay smooth. 
%   Breakpoints and engine tables are not changed.
%
%   The scope is a list of groups and/or single tables, each with its own
%   level. Later entries override earlier ones, so a general level can be
%   refined for a subset, for example
%
%     Create_Uncertainty(Data, {'all', 'control'}, [0.1, 0.3], seed)
%
%   gives +-10 % on all tables and +-30 % on the control derivatives.
%
%   Groups
%     'all'           All 22 aerodynamic tables
%     'none'          No uncertainty (nominal)
%     'longitudinal'  CL0 CLA CLDE CD0 CDA CM0 CMA CMDE CMQ
%     'lateral'       CYB CYDA CYDR CLLB CLLDA CLLDR CLLP CLLR
%                     CNB CNDA CNDR CNP CNR
%     'control'       CLDE CMDE CYDA CYDR CLLDA CLLDR CNDA CNDR
%     'static'        CLA CDA CMA CYB CLLB CNB
%     'damping'       CMQ CLLP CLLR CNP CNR
%     'lift'          CL0 CLA CLDE
%     'drag'          CD0 CDA
%     'sideforce'     CYB CYDA CYDR
%     'roll'          CLLB CLLDA CLLDR CLLP CLLR
%     'pitch'         CM0 CMA CMDE CMQ
%     'yaw'           CNB CNDA CNDR CNP CNR
%
%   Single tables, with or without '_tab'
%     Lift        CL0  CLA  CLDE
%     Drag        CD0  CDA
%     Side force  CYB  CYDA CYDR
%     Roll        CLLB CLLDA CLLDR CLLP CLLR
%     Pitch       CM0  CMA  CMDE CMQ
%     Yaw         CNB  CNDA CNDR CNP  CNR
%
%   One random number is drawn per table in a fixed order, whatever the
%   scope. The same seed therefore gives the same delta for a table in every
%   scope, so runs with different scopes can be compared directly.
%
%   Inputs
%     Data    struct Aero and engine data from Load_Aero_Data
%     scope   char or cell  Groups and/or table names (default 'all')
%     level   [1x1] or [Nx1] Maximum relative perturbation per scope entry,
%                            >= 0 (0.1 for +-10 %, default 0)
%     seed    [1x1]  Seed of the random stream (default 1), use a different
%                    seed per run, e.g. Sim.uncertainty.seed + k
%
%   Output
%     UNC   struct
%             Nominal    Unchanged Data
%             Perturbed  Data with perturbed aerodynamic tables, same layout,
%                        can be passed to the plant instead of Data
%             Delta      Applied perturbation per table, e.g. Delta.CLA_tab
%             Selected   Names of the perturbed tables
%             Scope      Scope entries and levels as given

%   Author:    Thomas Harald Müller
%   GitHub:    https://github.com/ThomasHaraldMueller
%   Copyright: (c) 2026 Thomas Harald Müller, Hamburg, Germany
%   Free to use and modify, provided this notice is kept and the author is credited.

%% Defaults and checks
if nargin < 2 || isempty(scope), scope = 'all'; end
if nargin < 3 || isempty(level), level = 0;     end
if nargin < 4 || isempty(seed),  seed  = 1;     end

scope = cellstr(scope);
level = level(:);
if isscalar(level)
    level = repmat(level, numel(scope), 1);
end

if numel(level) ~= numel(scope)
    error('Create_Uncertainty:level', ...
          'Give one level per scope entry or a single level (%d entries, %d levels).', ...
          numel(scope), numel(level));
end
if any(level < 0)
    error('Create_Uncertainty:level', 'Uncertainty levels must be non-negative.');
end

%% Tables, in a fixed order
Aero  = Data.AeroCoefficients;
names = sort(fieldnames(Aero));
names = names(endsWith(names, '_tab'));
n     = numel(names);

%% Random draw, one number per table independent of the scope
stream = RandStream('mt19937ar', 'Seed', seed);
r      = 2 * rand(stream, n, 1) - 1;               % Uniform in [-1, 1]

%% Level per table, later scope entries override earlier ones
groups    = coefficient_groups();
tab_level = zeros(n, 1);

for i = 1:numel(scope)
    sel            = resolve_scope(scope{i}, names, groups);
    tab_level(sel) = level(i);
end

%% Perturbation
UNC.Nominal = Data;
UNC.Delta   = struct();

for k = 1:n
    delta                = tab_level(k) * r(k);
    UNC.Delta.(names{k}) = delta;
    Aero.(names{k})      = Aero.(names{k}) * (1 + delta);
end

UNC.Perturbed                  = Data;
UNC.Perturbed.AeroCoefficients = Aero;
UNC.Selected                   = names(tab_level > 0);
UNC.Scope                      = struct('entry', scope(:), 'level', num2cell(level));

end

%% Local functions

function groups = coefficient_groups()
% COEFFICIENT_GROUPS  Table names (without '_tab') of each predefined group.
    groups.longitudinal = {'CL0','CLA','CLDE','CD0','CDA','CM0','CMA','CMDE','CMQ'};
    groups.lateral      = {'CYB','CYDA','CYDR','CLLB','CLLDA','CLLDR','CLLP','CLLR', ...
                           'CNB','CNDA','CNDR','CNP','CNR'};
    groups.control      = {'CLDE','CMDE','CYDA','CYDR','CLLDA','CLLDR','CNDA','CNDR'};
    groups.static       = {'CLA','CDA','CMA','CYB','CLLB','CNB'};
    groups.damping      = {'CMQ','CLLP','CLLR','CNP','CNR'};

    groups.lift         = {'CL0','CLA','CLDE'};
    groups.drag         = {'CD0','CDA'};
    groups.sideforce    = {'CYB','CYDA','CYDR'};
    groups.roll         = {'CLLB','CLLDA','CLLDR','CLLP','CLLR'};
    groups.pitch        = {'CM0','CMA','CMDE','CMQ'};
    groups.yaw          = {'CNB','CNDA','CNDR','CNP','CNR'};
end

function sel = resolve_scope(entry, names, groups)
% RESOLVE_SCOPE  Logical index of the tables covered by one scope entry.
    key = lower(strtrim(entry));

    if strcmp(key, 'all')
        sel = true(numel(names), 1);

    elseif strcmp(key, 'none')
        sel = false(numel(names), 1);

    elseif isfield(groups, key)
        tabs    = strcat(groups.(key), '_tab');
        missing = setdiff(tabs, names);
        if ~isempty(missing)
            error('Create_Uncertainty:group', ...
                  'Group ''%s'' refers to tables not in Data: %s', ...
                  key, strjoin(missing, ', '));
        end
        sel = ismember(names, tabs);

    else
        tab = [upper(regexprep(entry, '_tab$', '', 'ignorecase')), '_tab'];
        sel = strcmp(names, tab);
        if ~any(sel)
            error('Create_Uncertainty:scope', ...
                  ['Unknown scope entry ''%s''. Use a group (all, none, %s) ', ...
                   'or a table name.'], entry, strjoin(fieldnames(groups), ', '));
        end
    end
end