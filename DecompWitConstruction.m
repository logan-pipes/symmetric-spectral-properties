%%  DecompWitConstruction
%   Determines whether or not there exists a decomposable symmetric
%   entanglement witness with the maximal number (d(d-1)/2) of negative
%   eigenvalues.
%
%   This function has one required argument:
%     dimension: an integer >= 2
%
%   W = DecompWitConstruction(dimension) outputs
%   a decomposable symmetric entanglement witness with the maximal number
%   of negative eigenvalues, where the witness is of the form PV - (1/c)*P.
%   If no such matrix exists, it outputs NaN instead.
%
%   [W, c] = DecompWitConstruction(dimension) determines
%   whether there exists a decomposable symmetric entanglement
%   witness with the minimal number of negative eigenvalues that has
%   the form PV - (1/c)*P. It provides the witness and the value of c if
%   possible, or gives NaN for both if not possible.
%
%   [W, c, X, Y] = DecompWitConstruction(dimension)
%   makes the same calculation but provides the decomposition as well, with
%   X,Y such that W = PV*PartialTranspose(X)*PV + Y

%   requires: cvx (http://cvxr.com/cvx/), QETLAB (qetlab.com)
%   author: Logan Pipes (ldpipes@mta.ca)
%   last updated: August 20, 2021

%#ok<*EQEFF,*NODEF>

function [W, c, X, Y] = DecompWitConstruction(dimension)

if (dimension < 2 || ~(floor(dimension) == dimension))
    err = MException('DecompWitConstruction:invalidDimension', 'The dimension input to DecompWitConstruction must be an integer >= 2.');
    throw(err)
end

sym_rank = nchoosek(dimension,2); % Maximum d(d-1)/2 negative eigenvalues
spanning_mats = zeros(dimension,dimension,sym_rank);

counter = 0; % Generate #sym_rank symmetric linearly independent matrices
for j=1:dimension % whose forward diagonals all sum to zero
    for k=2:dimension-j+1
        counter = counter + 1;
        spanning_mats(1,j,counter) = 1;
        spanning_mats(j,1,counter) = 1;
        spanning_mats(k,k+j-1,counter) = -1;
        spanning_mats(k+j-1,k,counter) = -1;
    end
end

spanning_vecs = zeros(dimension^2,sym_rank);
for j=1:sym_rank
    % spanning_vecs is the span of the vectorizations of the previously
    % generated matrices
    spanning_vecs(:,j) = vec(spanning_mats(:,:,j));
end

% P is the orthogonal projection onto the columns of spanning_vecs
P = spanning_vecs * (spanning_vecs' * spanning_vecs)^-1 * spanning_vecs';

% PV is the projection onto the symmetric subspace
PV = SymmetricProjection(dimension);

cvx_begin sdp quiet
cvx_solver SDPT3
cvx_precision best
    variable X(dimension^2,dimension^2) hermitian semidefinite
    variable Y(dimension^2,dimension^2) hermitian semidefinite
    variable coef(1,1)
maximize coef
subject to
    W = PV - coef*P; % Has the maximum number of negative eigenvalues
    W == PV*PartialTranspose(X)*PV + Y; % Has to be decomposable
    X == PV*X*PV; % and X and Y have to be symmetric
    Y == PV*Y*PV;
cvx_end

if strcmp(cvx_status, 'Failed')
    W = NaN;
    c = NaN;
    X = NaN;
    Y = NaN;
else
    c = 1/coef;
end

end

