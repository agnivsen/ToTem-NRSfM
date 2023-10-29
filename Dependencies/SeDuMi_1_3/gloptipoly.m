function [output,sedumi] = gloptipoly(P,order,pars)
% GLOPTIPOLY  Global Optimization over Polynomials with SeDuMi
%
% [OUTPUT,SEDUMI] = GLOPTIPOLY(P,ORDER,PARS)
%
% Solve an LMI relaxation of a global optimization problem over
% polynomial and/or integer constraints
% 
% Input arguments:
%
% - P: cell array of polynomial constraints coefficients,
%   can take three formats:
%
%   - Each cell P{I}, I=1..M+1, is a structure with fields
%     - P{I}.c: coefficient matrix such that P{I}.c(J1,J2..JN) is the
%       coefficient of the monomial X(1)^(J1-1)*X(2)^(J2-1)*..*X(N)^(JN-1)
%       Notice that the coefficient matrix associated with a polynomial
%       with only one indeterminate must be a column vector
%     - P{I}.t: type of constraint, either '>=' (non-negative inequality),
%       '<=' (non-positive inequality), '==' (equality), 'min' (objective
%       function to minimize) or 'max' (objective function to maximize)
%     - P{I}.v: this optional field is a vector of integer constraints on
%       variables, where each entry P{I}.v(K) can take values:
%       - 0  meaning that variable X(K) is free (default)
%       - -1 meaning that variable X(K) must equal +1 or -1
%       - +1 meaning that variable X(K) must equal  0 or +1
%     This is a standard, explicit format for small-scale problems
%
%   - When a cell P{I} is a matrix and not a structure, then it is assumed
%     that P{I}.c = P{I}, P{I}.t = 'min' (objective function to minimize)
%     if I = 1 or P{I}.t = '>=' (non-negative inequality) if I > 1
%     This is an implicit format meant for notational simplicity
%
%   - When the row vector of maximum degrees [J1 J2 .. JN] is provided
%     as an additional field P{I}.s, then P{I}.c is assumed to be a
%     sparse column vector such that P{I}.c(SUB2IND(P{I}.s,J1,J2..JN) is the
%     coefficient of the monomial X(1)^(J1-1)*X(2)^(J2-1)*..*X(N)^(JN-1)
%     This is a format meant for medium-scale problems
%
%   - The trace of the moment matrix is minimized when no objective
%     function is specified in input structure P. To avoid this, enter a
%     constant objective function, for example P{I}.c = 0, P{i}.t = 'min'
%  
% - ORDER (optional): order of the LMI relaxation
%   It means that variables of degree up to 2*ORDER will be represented
%   An error message is issued when 2*ORDER is less than the largest
%   degree occuring in the problem
%   By default, 2*ORDER is equal to the largest problem degree
%
% - PARS (optional): parameter structure, including
%   - PARS.eps: desired accuracy in SeDuMi (default 1e-9)
%   - PARS.fid: 0 for no screen output (default 1)
%   - PARS.radius: feasibility radius in SeDuMi to prevent
%     unbounded problems. Primal SeDuMi variables are then constrained to
%     a Lorentz cone (default no feasibility radius)
%   - PARS.pert: scalar or vector of perturbation of the objective function
%     Useful to extract globally optimal solutions (default 0)
%   - PARS.scaling: scaling on decision variables. SeDuMi convergence is
%     greatly improved when optimal values of decision variables lie
%     around +1 (default no scaling)
%   - PARS.ranktol: threshold for evaluating ranks of moment matrices
%     when detecting global optimality and extracting solutions
%     (default 1e-3)
%   - PARS.pivotol: threshold for basis computation by Gaussian
%     elimination when extracting solutions (default 1e-6)
%   - PARS.testol: threshold for testing relaxed LMI solutions
%     (default 1e-6)
%  
% Output arguments:
%
% - OUTPUT: problem status, with solutions and objective at the
%   optimum when the global optimum was reached:
%   - OUTPUT.status: -1 if the relaxed LMI problem is infeasible or
%     could not be solved (see below the description of field SEDUMI.info)
%     - OUTPUT.obj is empty
%     - OUTPUT.sol  is empty
%   - OUTPUT.status: 0 if it is impossible to detect global optimality
%     - OUTPUT.obj is the optimum objective of the relaxed LMI problem
%     - OUTPUT.sol  is empty
%   - OUTPUT.status: 1 if the global optimum was reached, and
%     one or several isolated solutions were extracted
%     - OUTPUT.obj is the globally optimum objective
%     - OUTPUT.sol  is a cell array of vectors of globally optimal
%       decision variables
%   
% - SEDUMI: SeDuMi information variables
%   - SEDUMI.x, SEDUMI.y: optimal primal and dual SeDuMi variables
%   - SEDUMI.A, SEDUMI.b, SEDUMI.c: SeDuMi problem matrices
%   - SEDUMI.K: SeDuMi cone structure
%   - SEDUMI.info: SeDuMi information structure
%   - SEDUMI.M: cell array of moment matrices
%   - SEDUMI.pows: powers of the indeterminates
%   
%   The dual problem solved by SeDuMi is
%     max b'*y s.t. z = c-A'*y
%   where vector z belongs to the self dual-cone K made of K.f equality
%   constraints and K.s(1)..K.s(M+1) SDP constraints
%
%   Moment matrices of increasing sizes are stored in cell array
%   SEDUMI.M = {SEDUMI.M{1},..SEDUMI.M{ORDER}} where SEDUMI.M{1} is
%   the smallest upper-left (first order) moment matrix whose
%   first row or column contains relaxed polynomial variables
%   X(I) = SEDUMI.y(I), I=1..N
%
%   SEDUMI.pows is the matrix of variable powers such that
%   SEDUMI.y(I) = X(1)^SEDUMI.pows(I,1)*..*X(N)^SEDUMI.pows(I,N)
%   when the largest moment matrix SEDUMI.M{ORDER} has rank one
%
%   When SEDUMI.info.pinf = 1, then the SeDuMi primal problem is
%   infeasible and the LMI relaxation may be unbounded (see input
%   argument RADIUS)
%   When SEDUMI.info.dinf = 1, then the SeDuMi dual problem is
%   infeasible and the LMI relaxation may be infeasible as well
%   When SEDUMI.info.numerr > 0, then SeDuMi failed for numerical reasons
%
% The solver SeDuMi 1.05 must be installed

% Author: D. Henrion, January 24, 2001 (henrion@laas.fr)
% Last modified by D. Henrion, November 4, 2003
version = '2.2e';

addpath('D:\code\SEDUMI\sedumi-master\');

if exist('sedumi') ~= 2,
 error('SeDuMi is not properly installed');
end;

% ------------------------------------------------------------------------
% Check input arguments

if ~isa(P,'cell'),
  P = {P};
elseif min(size(P)) ~= 1,
  error('Invalid first input argument');
end;

% Number of constraints
nbconstr = length(P);

% Put P in standard format, cell array with components
% P{k}.c = matrix of polynomial coeffs, P{k}.t = type of constraint
typevar = [];
for k = 1:nbconstr,
 if ~isa(P{k}, 'struct'),
  if ~isa(P{k}, 'double'),
   error('Invalid structure of first input argument');
  else
   P{k}.c = P{k};
   if k == 1,
    P{k}.t = 'min'; % first P = objective function
   else
    P{k}.t = '>='; % geq constraint
   end;
  end;
 end;
 fn = fieldnames(P{k});
 c = []; t = []; s = [];
 for l = 1:length(fn),
  if fn{l}(1) == 'c',
   c = getfield(P{k}, fn{l});
  elseif fn{l}(1) == 't',
   t = getfield(P{k}, fn{l});
  elseif fn{l}(1) == 's',
   s = getfield(P{k}, fn{l});
  elseif fn{l}(1) == 'v',
   typevar = getfield(P{k}, fn{l});
  else
   error(['Invalid polynomial #' int2str(k)]);
  end;
 end;
 if isempty(c),
  error(['Invalid empty coefficient matrix in polynomial #' int2str(k)]);
 end;
 if isempty(t),
  P{k}.t = '>=';
 end;
 if ~isempty(s),
  if ~issparse(c) | (min(size(c)) > 1),
   error(['Coefficient must be a sparse vector in polynomial #' int2str(k)]);
  end;
 else
  P{k}.s = []; % not sparse
 end;
end;

% Number of variables, degrees of polynomials and size of coef matrices
nbvar = 0;
for k = 1:nbconstr,
 if isempty(P{k}.s), % not sparse
  if ndims(P{k}.c) == 2,
   if ~any(size(P{k}.c)),
    P{k}.s = 0; % empty
   elseif size(P{k}.c,2) == 1,
    P{k}.s = size(P{k}.c,1); % column vector = one variable
   else
    P{k}.s = size(P{k}.c); % row vector or 2D matrix = two variables
   end;
  else
   P{k}.s = size(P{k}.c); % more than 2D
  end;
 end;
 nbvar = max(nbvar, length(P{k}.s));
 degent = locate(find(P{k}.c),P{k}.s); % non-zero entries
 P{k}.d = 0;
 for l = 1:length(degent),
  P{k}.d = max(P{k}.d, sum(degent{l}-1));
 end;
end;
if nbvar < 1,
  error('No variables');
end;

% Parse constraint types and extract objective function
Pobj = []; indequ = []; indinequ = [];
for k = 1:nbconstr,
 switch P{k}.t,
  case {'>=', '=>'}, % inequality
   P{k}.t = 1;
   indinequ = [indinequ k];
  case {'<=', '=<'}, % inequality
   P{k}.t = 1;
   P{k}.c = -P{k}.c; % change sign
   indinequ = [indinequ k];
  case {'=='}, % equality
   P{k}.t = 0;
   indequ = [indequ k];
  case {'min', 'Min', 'MIN'},
   if isempty(Pobj),
    Pobj = P{k};
    signobj = 1;
   else
    error(['Objective defined twice in polynomial #' int2str(k)]);
   end;
  case {'max', 'Max', 'MAX'},
   if isempty(Pobj),
    Pobj = P{k};
    signobj = -1;
    Pobj.c = -Pobj.c; % change sign
   else    
    error(['Objective defined twice in polynomial #' int2str(k)]);
   end;
  otherwise
   error(['Invalid constraint type in polynomial #' int2str(k)]);
 end;
end;

% Sort constraints: first equalities, then moment matrix, then inequalities
Pmeas.c = 1; Pmeas.t = -1; Pmeas.s = []; Pmeas.d = 0; 
P = {P{:} Pmeas};
P = {P{[indequ nbconstr+1 indinequ]}};
nbconstr = length(P)-1;

% Maximum degree in objective and constraints
if ~isempty(Pobj),
 degpoly = Pobj.d;
else
 degpoly = 0;
end;
for k = 1:nbconstr+1,
 degpoly = max(degpoly, P{k}.d);
end;
if degpoly < 1,
  error('No variables');
end;
  
% Default input parameters
if nargin < 2, order = []; end;
if nargin < 3, pars = []; end;

if isempty(typevar),
 typevar = zeros(1,nbvar); % free variables
elseif ~isa(typevar,'double') | (length(typevar) ~= nbvar) | ...
      (min(size(typevar)) ~= 1),
  error('Invalid variable type vector');
end;
if ~isfield(pars,'radius'),
 pars.radius = -1; % no feasibility radius on decision variables
elseif ~isa(pars.radius,'double') | (length(pars.radius) ~= 1),
 error('Invalid feasibility radius');
end;
if ~isfield(pars,'fid'),
 pars.fid = 1; % trace on
end;
if ~isfield(pars,'eps'), 
 pars.eps = 1e-9; % SeDuMi accuracy
end;
if ~isfield(pars,'pert'),
 pars.pert = 0; % no pertubation of objective function
end;
if ~isfield(pars,'scaling'),
 pars.scaling = ones(1,nbvar); % no scaling of decision variables
end;
if ~isfield(pars,'ranktol'),
 pars.ranktol = 1e-3; % threshold for evaluating ranks of moment matrices
end;
if ~isfield(pars,'pivotol'),
 pars.pivotol = 1e-6; % threshold for basis computation via Gaussian elim
end;
if ~isfield(pars,'testol'),
 pars.testol = 1e-6; % threshold for testing LMI relaxed vector
end;

if isempty(order),
 % default LMI relaxation order =
 % half the maximum degree in constraints and objective
 order = ceil(degpoly/2);
elseif ~isa(order,'double'),
  error('Second input argument must be an integer');
else
 order = floor(order);
 if (length(order) > 1) | (order < 1),
  error('Second input argument must be a strictly positive integer');
 end;
end;

if pars.fid,
 disp(['GloptiPoly ' version ...
       ' - Global Optimization over Polynomials with SeDuMi']);
 if sum(abs(typevar)) == 0,
  disp(['Number of variables = ' int2str(nbvar)]);
 else
  disp(['Number of variables = ' int2str(nbvar) ...
        ' among which ' int2str(sum(abs(typevar))) ' are integers']);
 end;
 disp(['Number of constraints = ' int2str(nbconstr)]);
 disp(['Maximum polynomial degree = ' int2str(degpoly)]);
 disp(['Order of LMI relaxation = ' int2str(order)]);
 disp('Building LMI. Please wait..');
end;

% ********************************************************************
% Build LMI (decision) variables

if 2*order < degpoly,
 disp('Some variables are not represented in the LMI relaxation');
 error('Increase relaxation order');
end;

% Indices of LMI variables in base BASE=2*ORDER+1 such that
% DEC2BASE(INDICES,BASE) returns lexical indices of dual variables
indices = 0; base = 2*order+1;

if base^nbvar > 2^31-2,
 error(['Problem too large to be handled by the current version'...
       ' of Gloptipoly. Sorry']);
end;

degvarLMI = 0;
for k = 1:order,
 newindices = generate(nbvar,k,base);
 indices = [indices newindices];
 degvarLMI = [degvarLMI k*ones(1,length(newindices))];
end;

% Order of LMI = order of relaxation - degree of constraint / 2
for k = 1:nbconstr+1,
  P{k}.o = order - ceil(P{k}.d/2);
end;

% Matrix of variable products
basmat = zeros(length(indices));
for k = 1:length(indices),
  basmat(k,:) = indices(k)+indices;
end;

% Filter 0-1 and +/-1 variable constraints
basmat = filtermat(basmat,typevar,base);

% Remove redundant variables
retain = 1;
for k = 2:length(basmat),
 if ~any(basmat(1,1:k-1) == basmat(1,k)),
   retain = [retain k];
 end;
end;  
basmat = basmat(retain,retain);

% Number of LMI variables for each degree = size of LMI-1
% = index when defining LMI constraints
degvarLMI = degvarLMI(retain);
indLMI = zeros(1,order);
for k = 1:order,
 indLMI(k) = sum(degvarLMI == k);
end;

% Dimensions of LMIs
for k = 1:nbconstr+1,
 P{k}.dims = sum(indLMI(1:P{k}.o))+1;
end;

% Shift indices in SeDuMi matrices
% equality constraint = (dim+1)*dim/2 entries -> K.f = (dim+1)*dim/2
% inequality constraint = dim^2 entries -> K.s = dim
shift = zeros(1,nbconstr+1);
for k = 1:nbconstr+1,
 if P{k}.t,
  shift(k) = (P{k}.dims)^2;
 else
  shift(k) = (P{k}.dims+1)*(P{k}.dims)/2;
 end;
end;
shift = [0 cumsum(shift)];
sizeLMI = shift(end);

% Correspondance between lexical index and variable index
varLMIind = basmat(1,2); nbvarLMI = 1;
maxind = max(max(basmat));
lexind = sparse(maxind,1);
lexind(varLMIind) = nbvarLMI;
for row = 1:length(basmat),
 for col = row:length(basmat),
  newind = basmat(row,col);
  if (newind > 0) & ~any(varLMIind == newind), % new variable
   nbvarLMI = nbvarLMI + 1;
   lexind(newind) = nbvarLMI;
   varLMIind = [varLMIind newind];
  end;
 end;
end;

% Perturbation of the objective function
if ~isa(pars.pert,'double'),
 error('Invalid vector of perturbation of objective function');
elseif length(pars.pert) == 1, % scalar
 pars.pert = pars.pert*ones(1,nbvarLMI);
elseif length(pars.pert) == nbvar, % only original variables
 pars.pert = [pars.pert zeros(1,nbvarLMI-nbvar)];
elseif length(pars.pert) ~= nbvarLMI,
 error('Invalid vector of perturbation of objective function');
end;

% Scale problem if required
if ~all(pars.scaling == 1),
  if min(size(pars.scaling)) > 1,
    error('Incorrect scaling vector');
  end;
  if pars.fid,
    disp(['Scale problem. Scaling norm = ' num2str(norm(pars.scaling))]);
  end;
  for k = 0:length(P)
   if k > 0,
    coef = P{k}.c;
   elseif ~isempty(Pobj)
    coef = Pobj.c;
   else
    coef = [];
   end;
   if ~isempty(coef),
     if k > 0,
       sizeP = P{k}.s;
     else
       sizeP = Pobj.s;
     end;
     dimP = length(sizeP);
     if length(pars.scaling) == 1,
       pars.scaling = pars.scaling * ones(1,dimP);
     elseif length(pars.scaling) < dimP,
       error('Incorrect size of scaling vector');
     end;
     if ~issparse(coef), % non-sparse coefficient
       ind = cell(1,dimP);
       for i = 1:dimP,
	 ind{i} = 1:sizeP(i);
       end;
       for i = 1:dimP
	 for j = 2:sizeP(i)
	   newind = ind;
	   newind{i} = j;
	   coef(newind{:}) = coef(newind{:}) * pars.scaling(i)^(j-1);
	 end;
       end;
     else % sparse coefficient
       index = find(coef);
       power = locate(index,sizeP); % compute polynomial powers
       for i = 1:length(index),
	 coef(index(i)) = prod([coef(index(i)) pars.scaling.^(power{i}-1)]);
       end;
     end;
     if k > 0,
      P{k}.c = coef;
     else
      Pobj.c = coef;
     end;
    end;
  end;
end;

% ********************************************************************
% Build LMI matrices for SeDuMi
% in dual form: max b*y s.t. c-A*y in K
% as described in [Jean B. Lasserre. Global optimization with polynomials
% and the problem of moments. SIAM Journal on Optimization,
% Vol. 11, No. 3, pp. 796-817, 2001]

A = sparse(sizeLMI,nbvarLMI);
b = sparse(1,nbvarLMI);
c = sparse(sizeLMI,1);

for constr = 1:nbconstr+1,

 rowshift = shift(constr);

 if P{constr}.t < 0,

  % --------------
  % Moment matrix

  % locate occurrences of constant term in index matrix
  indmat = basmat;
  dim = length(indmat);
  for row = 1:dim,
   for col = row:dim,
    index = indmat(row,col);
    if index,
     var = lexind(index);
     % update LMI constraint matrix
     A(rowshift+row+(col-1)*dim, var) = -1;
     A(rowshift+col+(row-1)*dim, var) = -1;
    else
     % update right-hand side vector
     c(rowshift+row+(col-1)*dim) = 1;
     c(rowshift+col+(row-1)*dim) = 1;
    end;
   end;
  end;

 elseif P{constr}.t > 0,
   
  % ----------------------
  % Inequality constraints

  poly = P{constr}.c(:); % polynomial constraints into column vector
  nonzero = find(poly); poly = poly(nonzero); % extract non-zero components
  power = locate(nonzero,P{constr}.s); % compute polynomial powers

  refmat = basmat(1:P{constr}.dims,1:P{constr}.dims); % index ref matrix
  for term = 1:length(nonzero), % for each term in the polynomial
   % compute index matrix
   if max(power{term}) > base,
    disp('Some variables are not represented in the constraint matrix');
    error('Increase relaxation order');
   end;
   indvar = 0;
   for k = 1:length(power{term}),
    indvar = indvar+(power{term}(k)-1)*base^(k-1);
   end;
   % indvar = base2dec(fliplr(sprintf('%d',power{term}-1)),base);
   % shift reference matrix and filter 0-1 and +/-1 constraints
   indmat = filtermat(refmat+indvar, typevar, base);
   dim = length(indmat);
   for row = 1:dim,
    for col = row:dim,
     index = indmat(row,col);
     if index,
      if index <= maxind,
       var = lexind(index);
       if ~var,
        disp('Some variables are not represented in the constraint matrix');
        error('Increase relaxation order');
       end;
      else
       disp('Some variables are not represented in the constraint matrix');
       error('Increase relaxation order');
      end;
      % update constraint matrix
      A(rowshift+row+(col-1)*dim, var) = -poly(term);
      A(rowshift+col+(row-1)*dim, var) = -poly(term);
     else
      % update right-hand side vector
      c(rowshift+row+(col-1)*dim) = poly(term);
      c(rowshift+col+(row-1)*dim) = poly(term);     
     end;
    end;
   end;
  end;
 
 else,
   
  % --------------------
  % Equality constraints

  poly = P{constr}.c(:); % polynomial constraints into column vector
  nonzero = find(poly); poly = poly(nonzero); % extract non-zero components
  power = locate(nonzero,P{constr}.s); % compute polynomial powers

  refmat = basmat(1:P{constr}.dims,1:P{constr}.dims); % index ref matrix
  for term = 1:length(nonzero), % for each term in the polynomial
   % compute index matrix
   if max(power{term}) > base,
    disp('Some variables are not represented in the constraint matrix');
    error('Increase relaxation order');
   end;
   indvar = 0;
   for k = 1:length(power{term}),
    indvar = indvar+(power{term}(k)-1)*base^(k-1);
   end;
   % indvar = base2dec(fliplr(sprintf('%d',power{term}-1)),base);
   % shift reference matrix and filter 0-1 and +/-1 constraints
   indmat = filtermat(refmat+indvar, typevar, base);
   dim = length(indmat);
   indrow = rowshift;
   for row = 1:dim,
    for col = row:dim,
     index = indmat(row,col);
     indrow = indrow+1;
     if index,
      if index <= maxind,
       var = lexind(index);
       if ~var,
        disp('Some variables are not represented in the constraint matrix');
        error('Increase relaxation order');
       end;
      else
       disp('Some variables are not represented in the constraint matrix');
       error('Increase relaxation order');
      end;
      % update constraint matrix
      A(indrow, var) = -poly(term);
     else
      % update right-hand side vector
      c(indrow) = poly(term);
     end;
    end;
   end;
  end;
   
 end; % if

end;

% ---------
% Objective

if ~isempty(Pobj),
  
 poly = Pobj.c(:); % polynomial objective into column vector
 nonzero = find(poly); poly = poly(nonzero); % extract non-zero components
 power = locate(nonzero,Pobj.s); % compute polynomial powers

 constant = 0; % constant term
 for term = 1:length(nonzero), % for each term in the polynomial
  % compute index variable
  indvar = 0;
  for k = 1:length(power{term}),
   indvar = indvar+(power{term}(k)-1)*base^(k-1);
  end;
  % indvar = base2dec(fliplr(sprintf('%d',power{term}-1)),base);
  % filtered index variable
  index = filtermat(indvar, typevar, base);
  if index,
   if index <= maxind,
    var = lexind(index);
    if ~var,
     disp('Some variables are not represented in the objective function');
     error('Increase relaxation order');
    end;
   else
    disp('Some variables are not represented in the objective function');
    error('Increase relaxation order');
   end;  % update objective function vector
   b(var) = -poly(term);
  else
   % constant term
   constant = constant-poly(term);
  end;
 end;

else
  
 % If there is no objective function to be optimized, then
 % minimize the trace of the moment matrix
 for k = 1:length(basmat),
  index = basmat(k,k);
  if index,
   var = lexind(index);
   b(var) = -1;
  end;
 end;

end;

% Perturbation of objective function
if norm(pars.pert) > 0,
 b = b + pars.pert;
end;

% SeDuMi structure of LMI
K.f = 0; K.s = [];
for k = 1:nbconstr+1,
 if P{k}.t,
  K.s = [K.s P{k}.dims]; % inequality
 else
  K.f = K.f+(P{k}.dims+1)*(P{k}.dims)/2; % equality
 end;
end;

% Add quadratic cone constraint on vector of decision variables
% to prevent unbounded solutions
if pars.radius > 0, % set feasibility radius
 nbeq = sum(K.f);
 A = [A(1:nbeq,:); sparse(1,nbvarLMI); speye(nbvarLMI); A(nbeq+1:end,:)];
 c = [c(1:nbeq,:); pars.radius; sparse(nbvarLMI,1); c(nbeq+1:end,:)];
 K.q = nbvarLMI+1;
elseif pars.radius == 0, % unrestricted feasibility radius
 nbeq = sum(K.f);
 A = [A(1:nbeq,:); sparse(1,nbvarLMI); speye(nbvarLMI); A(nbeq+1:end,:)];
 c = [c(1:nbeq,:); 0; sparse(nbvarLMI,1); c(nbeq+1:end,:)];
 K.q = nbvarLMI+1; 
 nbvarLMI = nbvarLMI+1;
 A(nbeq+1,nbvarLMI) = 1;
 b(nbvarLMI) = 0;
end;

% Display information
if pars.fid,
 disp(['Number of LMI decision variables = ' int2str(nbvarLMI)]);
 disp(['Size of LMI constraints = ' int2str(sizeLMI)]);
 disp(['Sparsity of LMI constraints = ' ...
      num2str(100*nnz(A)/prod(size(A))) '% of non-zero entries']);
 disp(['Norm of perturbation of objective function = ' ...
      num2str(norm(pars.pert))]);
 disp(['Numerical accuracy for SeDuMi = ' num2str(pars.eps)]);
 if pars.radius > 0,
  disp(['Feasibility radius = ' num2str(pars.radius)]);
 elseif pars.radius == 0,
  disp('Unrestricted feasibility radius');
 else
  disp('No feasibility radius');   
 end;
 disp('Solving LMI problem with SeDuMi..');
end;

% ********************************************************************
% Call SeDuMi to solve LMI

save data A b c K
[x,y,info] = sedumi(A,b,c,K,pars);

% Diagnostic
status = 0;

if pars.fid,
  disp(['CPU time = ' num2str(info.cpusec ) ' sec']);
end;

if info.numerr == 2,
 if pars.fid,
  if pars.radius < 0,
   disp('SeDuMi dual problem may be unbounded');
   disp('Try to enforce feasibility radius');
  end;
  disp('Numerical problems');
 end;
 obj = [];
elseif info.pinf | info.dinf,
 if pars.fid,
  if info.pinf,
   if pars.radius < 0,
    disp('SeDuMi primal problem is infeasible');
    disp('SeDuMi dual problem may be unbounded');
    disp('Try to enforce feasibility radius');
   else
    disp('SeDuMi primal problem is infeasible');
   end;
  elseif info.dinf,
   disp('SeDuMi dual problem is infeasible');
   if pars.radius > 0,
    disp('Try to enlarge feasibility radius');
   else
    disp('Try to scale problem variables');
    disp('Original optimization problem may be infeasible');
   end;
  end;
 end;
 obj = [];
else,
 status = 1;
end;

% ------------------------------------------------------------------------
% Retrieve variables

% Norm of vector of decision variables
normy = norm(y);
if normy > 1e9,
 status = 0;
  if pars.fid,
   disp(['Norm of SeDuMi dual vector = ' num2str(normy) ' > 1e9']);
   disp('Warning ! SeDuMi dual may be unbounded');
   disp('Try to enforce feasibility radius or to scale the problem')
  end;
end;

% SeDuMi output vector
sedumi.x = x; sedumi.y = y; sedumi.info = info;
sedumi.A = A'; sedumi.b = b'; sedumi.c = c; sedumi.K = K;
sedumi.pars = pars;

% LMI sizes
% sedumi.dims = zeros(1,nbconstr+1);
% for i = 1:nbconstr+1,
%  if P{i}.t, % inequality
%  sedumi.dims(i) = -P{i}.dims;
% else % equality
%  sedumi.dims(i) = P{i}.dims;
% end;
% end;

% Variable types
% sedumi.v = typevar;

% Identification vector
ind = varLMIind';
sedumi.pows = zeros(size(ind,1),nbvar);
for k = 1:nbvar,
 sedumi.pows(:,k) = rem(ind,base);
 ind = (ind-sedumi.pows(:,k))/base;
end;

% Extract moment matrices of increasing sizes
z = sedumi.c-sedumi.A'*sedumi.y;
m = sedumi.K.s(1); MM = zeros(m);
MM(:) = z(sedumi.K.f+[1:m^2]);
sp = sum(sedumi.pows,2);
sedumi.M = {}; d = 0; k = 0;
while d < m,
 k = k+1; d = 1+sum(sp<=k);
 sedumi.M{k} = MM(1:d,1:d);
end;

if status % if LMI was solved properly
  
  % Compute objective function
  if ~isempty(Pobj),
    obj = signobj*(-constant-b*y);
    if pars.fid,
      disp(['LMI objective = ' num2str(obj)]);
    end;
  else
    obj = [];
    if pars.fid,
      disp('No objective function. Trace of moment matrix was minimized');
    end;
  end;

  % Check relaxed LMI vector: if it reaches the LMI objective
  % and it is feasible, then the global optimum was reached
  globopt = 0;
  reach = 0; feas = 0;
  if pars.fid,
    disp(['Checking relaxed LMI vector with threshold = ' ...
	 num2str(pars.testol)]);
  end;
  x = {sedumi.y(1:nbvar)};
  
  % objective function
  if ~isempty(Pobj),
    val = evaluate(Pobj, x{1});
    reach = abs(val-obj) < pars.testol * abs(val);
    if pars.fid,
      if reach,
	disp('Relaxed vector reaches LMI objective');
      else
	disp(['Relaxed vector reaches objective ' ...
	      num2str(val)]);
      end;
    end;
  else
    reach = 1; % no objective function
  end;
  feas = 1; k = 0;
  while feas & (k < nbconstr+1),
    k = k+1;
    if P{k}.t > 0 % inequality
      val = evaluate(P{k}, x{1});
      feas = val > -pars.testol;
    elseif P{k}.t == 0 % equality
      val = evaluate(P{k}, x{1});
      feas = abs(val) < pars.testol;
    end;
  end;
  if pars.fid,
    if feas,
      disp('Relaxed vector is feasible');
    else
      disp('Relaxed vector is not feasible');
    end;
  end;
  globopt = reach & feas;

  if globopt,
    
    % **************************
    % Relaxed vector is solution
    
    if pars.fid,
      disp('Global optimum was reached');
    end;

  else
    
    % *****************************************
    % Compute SVD of moment matrices
    % to detect global optimality

    % half of maximum degree - rank shift to detect global optimum
    rankshift = 1;
    for k = 1:nbconstr+1,
      rankshift = max(rankshift, ceil(P{k}.d/2));
    end;
    
    if pars.fid,
      disp(['Detecting global optimality (rank shift = ' ...
	    int2str(rankshift) ')..']);
      disp(['Relative threshold for rank evaluation  = ' ...
	    num2str(pars.ranktol)]);
    end;
    
    kmax = length(sedumi.M);
    U = cell(1,kmax); S = cell(1,kmax);
    rankM = zeros(1,kmax);
    
    oldrank = 1;
    rankdiff = 0;
    
    k = 0;
    while (k < kmax) & (globopt == 0),
      k = k + 1;
      [Uk,Sk,Vk] = svd(sedumi.M{k});
      U{k} = Uk; S{k} = Sk;
      S{k} = diag(S{k}); n = length(S{k});
      drop = find(S{k}(2:n) ./ S{k}(1:n-1) < pars.ranktol);
      if ~isempty(drop), rankM(k) = drop(1); else rankM(k) = n; end;
      if pars.fid,
	fprintf(['Moment matrix of order %2d has size %2d ' ...
		 'and rank %2d\n'], k, size(U{k},1),rankM(k));
	
      end;
      if rankM(k) <= oldrank,
	rankdiff = rankdiff + 1;
      else
	rankdiff = 0;
      end;
      if (rankdiff >= rankshift) & (globopt == 0),
	globopt = k; % global optimum was reached
	if pars.fid,
	  disp('Rank condition ensures global optimality');
	end;
      end;
      oldrank = rankM(k);
    end;
    
    % solutions
    x = {};
    
    % ********************************************************************
    % Extract solutions
    %
    % Algorithm:
    % 1. extract a Cholesky factor U of moment matrix SEDUMI.M{K} such that
    %    U*U' = SEDUMI.M{K}
    % 2. reduce U to column echelon form via Gaussian elimination with
    %    column pivoting, identify monomial basis B(1)..B(P) in U (pivots)
    % 3. for each monomial X(I), I=1..P extract the coefficient
    %    matrix N(I) of monomials X(I)*B(1)..X(I)*B(P) in basis B(1)..B(P)
    % 4. compute common roots of multiplication matrices N(1)..N(P) and store
    %    them in output cell array X
    
    if pars.fid,
      disp('Extracting solutions..');
      disp(['Relative threshold for basis detection  = ' ...
	    num2str(pars.pivotol)]);
    end;

    % Step 1: extract Cholesky factor of moment matrix
    U{k} = U{k}(:,1:rankM(k))*diag(sqrt(S{k}(1:rankM(k))));
    
    % Step 2: reduce Cholesky factor to column echelon form
    % and identify monomial basis
    [U{k},basis] = cef(U{k},pars.pivotol);

    % Step 3: extract multiplication matrix for each variable
    N = cell(1,nbvar);
    mon = [zeros(1,nbvar); sedumi.pows(1:n-1,:)]; % powers of monomials
    nmon = size(mon,1); % total number of monomials

    i = 0; fail = 0;
    while ~fail & (i < nbvar),
      i = i + 1; j = 0; M = zeros(rankM(k));
      monvar = zeros(1,nbvar); monvar(i) = 1;
      while ~fail & (j < rankM(k)),
	j = j + 1;
	% find monomial X(I)*B(J)
	% filter 0-1 and +/-1 constraints (cf. GLOPTIPOLY)
	newmon = filtermat(monvar+mon(basis(j),:),typevar,2*k+1);
	row = find(all(mon==ones(nmon,1)*newmon,2));
	if isempty(row),
	  % not represented in the moment matrix, so increase order
	  fail = 1;
	else
	  % build coefficient matrix
	  M(j,:) = U{k}(row,:);
	end;
      end;
      if ~fail,
	N{i} = M;
      end;
    end;
    
    if ~fail,

      % Step 4: compute common roots with the algorithm described in
      % [Corless, Gianni, Trager. A reordered Schur factorization method
      % for zero-dimensional polynomial systems with multiple roots.
      % Proc. ISSAC (Maui), W. Kuechlin (Ed) pp. 133-140, 1997]
	
      % random combinations of multiplication matrices
      coef = rand(nbvar,1); coef = coef / sum(coef);
      M = zeros(rankM(k));
      for i = 1:nbvar,
	M = M + coef(i)*N{i};
      end;
	
      % ordered Schur decomposition of M      
      [Q,T] = orderschur(M);
	
      % retrieve optimal vectors
      % it is assumed than there is no multiple root
      x = cell(1,rankM(k));
      maxerr = -Inf;
      for i = 1:rankM(k),
	x{i} = zeros(nbvar,1);
	for j = 1:nbvar,
	  x{i}(j) = Q(:,i)'*N{j}*Q(:,i);
	end;
	u = prod((ones(nmon,1)*x{i}').^mon, 2); % evaluate monomials
	err = max(abs(U{k}*u(basis)-u)) / max(abs(u));
	if err > maxerr,
	  maxerr = err;
	end;
      end;
      
      if pars.fid,
	disp(['Maximum relative error = ' num2str(maxerr)]);
      end;

      if 0
      % check if extracted solutions reach objective function
      reach = 1;
      if ~isempty(Pobj)
       for i = 1:rankM(k)  
        val = evaluate(Pobj, x{i});
        reach = reach & (abs(val-obj) < pars.testol * abs(val));
       end
       if pars.fid,
        if reach,
         disp('Extracted vectors reaches LMI objective')
        else
	 disp('Some extracted vectors do not reach LMI objective')
        end;
       end
      end
      % check if extracted solutions are feasible      
      feas = 1;
      for i = 1:rankM(k)
       k = 0;
       while feas & (k < nbconstr+1),
        k = k+1;
        if P{k}.t > 0 % inequality
         val = evaluate(P{k}, x{i});
         feas = val > -pars.testol;
        elseif P{k}.t == 0 % equality
         val = evaluate(P{k}, x{i});
         feas = abs(val) < pars.testol;
        end;
       end;
      end
      if pars.fid,
       if feas,
        disp('Extraxed vectors are feasible');
       else
        disp('Some extracted vectors are not feasible');
       end;
      end;

      globopt = reach & feas;
      stop = 1;
     end
     
    else
      
      if pars.fid,
	fprintf('Incomplete basis\n');
      end;
      
    end;
    
  end;
  
  % number of solutions
  nbsol = length(x);
  
  % scale back variables
  if nbsol > 0
    if ~all(pars.scaling == 1),
      if pars.fid,
	disp('Scale back decision variables');
      end;
      for k = 1:nbsol,
	for i = 1:nbvar,
	  x{k}(i) = x{k}(i) * pars.scaling(i);
	end;
      end;
    end;
  end;

end; % if status > 0

% ***********************
% Assign output argument

if status == 0 % LMI problem has no solution or could not be solved
  output.status = -1; 
  output.obj = [];
  output.sol = {};
elseif globopt == 0 % impossible to dectect global optimality
  output.status = 0;
  output.obj = obj;
  output.sol = {};
  if pars.fid,
    disp('Impossible to detect global optimality');
    if ~isempty(Pobj),
      if signobj,
        disp('LMI objective is a lower bound on the global minimum');
      else
	disp('LMI objective is an upper bound on the global minimum');
      end;
    end;
  end;
else % global optimum was reached and solutions were extracted
  output.status = 1;
  output.obj = obj;
  output.sol = x;
  if pars.fid,
   if nbsol > 1,
    disp([int2str(nbsol) ' solutions extracted']);
   else
    disp('One solution extracted');
   end;
  end;
end;

return;

% ********************************************************************
% Utilities

function mat = filtermat(mat,typevar,base)
% Filter 0-1 and +/-1 constraints in matrix MAT
% according to type variable vector TYPEVAR and base BASE
if any(typevar),
 nbvar = length(typevar);
 % Decomposition of MAT in base BASE
 mask(:,:,nbvar) = mat;
 for k = nbvar:-1:2,
   remainder = rem(mask(:,:,k), base);
   mask(:,:,k-1) = (mask(:,:,k)-remainder)/base;
   mask(:,:,k) = remainder;
 end;
 % Recomposition with filtering
 mat = zeros(size(mat));
 for k = 1:nbvar,
  mat = mat * base;
  switch typevar(k),
   case 1, % 0-1 constraint : x^2 -> x
    mat = mat + (mask(:,:,k) ~= 0); 
   case -1, % +/-1 constraint: x^2 -> 0
    mat = mat + rem(mask(:,:,k), 2);
   case 0, % no constraint
    mat = mat + mask(:,:,k);
   otherwise
    error('Invalid entry in variable type vector');
  end;
 end;
end;

function multi = locate(indices,size)
% Transform linear indices INDICES into a cell array of
% multidimensional indices in a matrix of size SIZE
% For example, LOCATE(5,[2 3]) returns [1 3], i.e. 5th element in
% 2x3 matrix is located in first row and third column 
if any(size),
 multi = cell(length(indices),1);
 for k1 = 1:length(indices),
  index = indices(k1);
  for k2 = 1:length(size),
   multi{k1}(k2) =  1+rem(index-1,size(k2));
   index = 1+floor((index-1)/size(k2));
  end;
 end;
else
 multi = {1};
end;

function t = generate(digits,sum,base)
% Generate all the indices with DIGITS digits summing to SUM in base BASE.
% For example, GENERATE(3,2,3) returns the vector [2 4 10 6 12 18] which is
% [002 011 101 020 110 200] in base 3. Note that the most significant
% digit is located at the left.  
if digits < 1,
  t = [];
elseif digits < 2,
  t = sum;
else,
  t = zeros(1,nchoosek(sum+digits-1,digits-1));
  j = 1;
  for i = sum:-1:0,
    s = generate(digits-1,sum-i,base); r = length(s); % recursive call
    t(j:j+r-1) = i*ones(1,r)+s*base; j = j+r;
  end;
end;

function [A,basis] = cef(A,tol)
% The instruction
%
%  [E,BASIS] = CEF(A)
%
% computes a column echelon form E of matrix A
% and returns basis row indices in vector BASIS
%
% The relative threshold for column pivoting can be specified
% as an additional input argument
%
% The reduction is performed by Gaussian elimination with
% column pivoting, based on Matlab's RREF routine
  
[n,m] = size(A);

% Loop over the entire matrix.
i = 1; j = 1; basis = [];
while (i <= m) & (j <= n)
   % Find value and index of largest element in the remainder of row j
   [p,k] = max(abs(A(j,i:m))); k = k+i-1;
   if (p <= tol)
      % The row is negligible, zero it out.
      A(j,i:m) = zeros(1,m-i+1,1);
      j = j + 1;
   else
      % Remember row index
      basis = [basis j];
      % Swap i-th and k-th columns
      A(j:n,[i k]) = A(j:n,[k i]);
      % Find a non-negligible pivot element in the column
      found = 0;
      while ~found,
	if abs(A(j,i)) < tol*max(abs(A(:,i)))
	  j = j + 1;
	  found = (j == n);
	else
	  found = 1;
	end;
      end;
      if j <= n,
        % Divide the pivot column by the pivot element
	A(j:n,i) = A(j:n,i)/A(j,i);
	% Subtract multiples of the pivot column from all the other columns
	for k = [1:i-1 i+1:m]
	  A(j:n,k) = A(j:n,k) - A(j,k)*A(j:n,i);
	end
	i = i + 1;
	j = j + 1;
      end
   end
end

function P = evaluate(P,x)
% The instruction VAL = EVALUATE(P,X) evaluates an objective function,
% inequality or equality constraint P (in GloptiPoly's format) at a
% vector X

if issparse(P.c) % sparse coefficient matrix

 nd = length(P.s);
 d = P.s;
 val = 0;
 nzndx = find(P.c);
 for j = 1:length(nzndx),
  k = [1 cumprod(d(1:nd-1))];
  ndx = nzndx(j) - 1;
  pow = zeros(1,nd);
  for i = nd:-1:1,
   pow(i) = floor(ndx/k(i));
   ndx = rem(ndx,k(i));
  end;
  coef = 1;
  for i = 1:nd,
   coef = coef*x(i)^pow(i);
  end;
  val = val + P.c(nzndx(j))*coef;
 end;
 P = val;
 
else % non-sparse matrix, applies Horner's scheme

 P = P.c;
 nd = ndims(P);
 if (nd == 2) & (size(P,2) == 1), nd = 1; end;
 d = size(P);
 for n = 1:nd, ind{n} = 1:d(n); end;
 for n = nd:-1:1,
  indP = {ind{1:n-1} d(n)};
  newP = P(indP{:});
  for k = d(n)-1:-1:1,
   indP = {ind{1:n-1} k};
   newP = P(indP{:}) + x(n)*newP;
  end;
  P = newP;
 end;

end;

function [U,T] = orderschur(X)
% [U,T] = ORDERSCHUR(X) computes the real Schur decomposition X = U*T*U'
% of a matrix X with real eigenvalues sorted in increasing order along
% the diagonal of T
%
% Algorithm: perform unordered Schur decomposition with Matlab's SCHUR
% function, then order the eigenvalues with Givens rotations, see
% [Golub, Van Loan. Matrix computations. 1996]

[U,T] = schur(X);
U = real(U); T = real(T);
n = size(X,1);

order = 0;
while ~order, % while the order is not correct
 order = 1;
 for k = 1:n-1,
  if T(k,k) - T(k+1,k+1) > 0,
   order = 0; % diagonal elements to swap
   % Givens rotation
   [c,s] = givens(T(k,k+1),T(k+1,k+1)-T(k,k));
   T(k:k+1,k:n) = [c s;-s c]'*T(k:k+1,k:n);
   T(1:k+1,k:k+1) = T(1:k+1,k:k+1)*[c s;-s c];
   U(1:n,k:k+1) = U(1:n,k:k+1)*[c s;-s c];
  end;
 end;
end; % while

function [c,s] = givens(a,b)
% Givens rotation for ordered Schur decomposition
if b == 0
  c = 1; s = 0;
else
  if abs(b) > abs(a)
    t = -a/b; s = 1/sqrt(1+t^2); c = s*t;
  else
    t = -b/a; c = 1/sqrt(1+t^2); s = c*t;
  end
end






