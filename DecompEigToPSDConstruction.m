﻿%%  DecompEigToPSDConstruction
%   Uses semidefinite programming to implement Theorem 4 of the paper, and
%   thus shows that certain sets of numbers cannot be spectra of
%   decomposable symmetric entanglement witnesses.
%   
%   This function has one required argument:
%     mu: a row or column vector of d(d+1)/2 real eigenvalues, for some
%         integer 1 <= d <= 5.
%
%   pos = DecompEigToPSDConstruction(mu) outputs a boolean implementing
%   Theorem 4 from the paper: pos = 1 if there exist matrices Y such that
%   the inequalities all hold, and pos = 0 if no such matrices Y exist (so
%   there does not exist a decomposable symmetric entanglement witness with
%   eigenvalues mu).
%
%   [pos, Y] = DecompEigToPSDConstruction(mu) performs the same
%   calculation, but also outputs a cell array of k d*d matrices, for some
%   integer  1 <= k <= (d(d+1)/2)! such that if pos = 1 then they are the Y
%   matrices from Theorem 4, or NaN if no such set of matrices exist
%   (i.e., if pos = 0).

%   requires: cvx (http://cvxr.com/cvx/)
%   author: Logan Pipes (ldpipes@mta.ca)
%   last updated: August 20, 2021

%#ok<*AGROW,*NOPRT,*EQEFF>

function [pos, Y] = DecompEigToPSDConstruction(mu)

% mu is a row vector containing the eigenvalues in non-increasing order
mu = sort(mu(:)', 'descend'); %#ok<UDIM>

[h,len] = size(mu);
d = floor(sqrt(2*len));

if ~(h == 1) || ~(d*(d+1)/2 == len) || ~isreal(mu)
    err = MException('DecompEigToPSDConstruction:invalidDimension', 'The mu input to DecompEigToPSDConstruction must be a real vector of d(d+1)/2 eigenvalues, for some integer 1 <= d <= 5.');
    throw(err)
end

% We already know that the len! permutations simplify to just a few in the
% 2 <= d <= 5 cases, and have hard coded those in as follows
if d == 2
    matricizations = [3 1 2];

elseif d == 3
    matricizations = [6 2 1 3 4 5];

elseif d == 4
    matricizations = [10 8 6 1 7 5 2 4 3 9;
                      10 9 6 1 8 5 2 4 3 7;
                      10 9 7 1 8 6 2 4 3 5;
                      10 9 3 1 7 4 2 5 6 8;
                      10 9 3 1 8 4 2 5 6 7;
                      10 3 2 1 4 5 6 7 8 9;
                      10 9 2 1 8 4 3 5 6 7;
                      10 9 2 1 7 4 3 5 6 8;
                      10 8 7 1 6 5 2 4 3 9;
                      10 9 2 1 6 4 3 5 7 8;
                      10 6 2 1 5 4 3 7 8 9;
                      10 7 2 1 5 4 3 6 8 9;
                      10 3 2 1 4 5 7 6 8 9;
                      10 9 8 1 7 6 2 5 3 4;
                      10 8 2 1 6 4 3 5 7 9;
                      10 9 6 1 7 5 2 4 3 8;
                      10 8 2 1 5 4 3 6 7 9;
                      10 9 8 1 6 5 2 4 3 7;
                      10 8 3 1 7 4 2 5 6 9;
                      10 9 2 1 5 4 3 6 7 8;
                      10 8 3 1 6 4 2 5 7 9;
                      10 9 7 1 8 6 2 5 3 4;
                      10 9 7 1 6 5 2 4 3 8;
                      10 9 8 1 7 5 2 4 3 6;
                      10 9 8 1 7 6 2 4 3 5;
                      10 9 7 1 8 5 2 4 3 6];

elseif d == 5
    matricizations = [15,14,11,3,1,13,10,5,2,8,6,4,7,9,12;
                      15,14,13,9,1,12,11,8,2,10,7,3,6,4,5;
                      15,13,12,4,1,11,10,5,2,9,6,3,7,8,14;
                      15,13,4,2,1,11,6,5,3,7,8,9,10,12,14;
                      15,14,5,3,1,13,6,4,2,7,8,10,9,11,12;
                      15,14,12,8,1,13,11,7,2,10,6,3,5,4,9;
                      15,14,11,4,1,12,10,5,2,9,6,3,7,8,13;
                      15,4,3,2,1,5,6,7,9,8,10,11,12,13,14;
                      15,12,4,2,1,9,6,5,3,7,8,10,11,13,14;
                      15,13,11,3,1,10,8,5,2,7,6,4,9,12,14;
                      15,13,5,2,1,10,6,4,3,7,8,9,11,12,14;
                      15,14,10,3,1,13,9,4,2,7,6,5,8,11,12;
                      15,14,12,3,1,13,10,4,2,7,6,5,8,9,11;
                      15,14,4,3,1,12,6,5,2,7,8,10,9,11,13;
                      15,14,11,2,1,13,8,4,3,7,6,5,9,10,12;
                      15,9,3,2,1,7,6,5,4,8,10,11,12,13,14;
                      15,14,11,9,1,12,10,8,2,7,6,3,5,4,13;
                      15,14,5,3,1,13,6,4,2,7,8,9,10,11,12;
                      15,14,13,2,1,12,11,5,3,10,6,4,7,8,9;
                      15,14,13,2,1,9,8,5,3,7,6,4,10,11,12;
                      15,14,11,3,1,13,10,4,2,8,6,5,7,9,12;
                      15,14,9,2,1,11,8,4,3,7,6,5,10,12,13;
                      15,14,12,4,1,13,10,5,2,9,6,3,7,8,11;
                      15,11,5,2,1,10,6,4,3,7,8,9,12,13,14;
                      15,12,11,2,1,9,8,5,3,7,6,4,10,13,14;
                      15,4,3,2,1,5,6,7,11,8,9,12,10,13,14;
                      15,14,11,3,1,13,10,4,2,7,6,5,8,9,12;
                      15,14,12,10,1,13,11,9,2,8,7,3,6,4,5;
                      15,12,5,2,1,10,6,4,3,7,8,9,11,13,14;
                      15,14,13,9,1,12,10,7,2,8,6,3,5,4,11;
                      15,14,12,10,1,13,11,9,2,8,7,3,5,4,6;
                      15,14,11,4,1,13,10,5,2,9,6,3,7,8,12;
                      15,13,12,8,1,11,10,7,2,9,6,3,5,4,14;
                      15,14,12,10,1,13,11,8,2,9,7,3,6,4,5;
                      15,14,12,2,1,13,8,4,3,7,6,5,9,10,11;
                      15,14,11,3,1,12,9,4,2,7,6,5,8,10,13;
                      15,8,3,2,1,7,6,5,4,9,10,11,12,13,14;
                      15,14,12,3,1,13,11,4,2,8,6,5,7,9,10;
                      15,14,5,3,1,12,6,4,2,7,8,9,10,11,13;
                      15,13,9,3,1,12,8,4,2,7,6,5,10,11,14;
                      15,8,3,2,1,7,6,5,4,9,10,12,11,13,14;
                      15,14,13,3,1,12,11,5,2,9,6,4,7,8,10;
                      15,13,5,2,1,11,6,4,3,7,8,9,10,12,14;
                      15,11,9,2,1,10,8,4,3,7,6,5,12,13,14;
                      15,14,13,4,1,11,10,5,2,9,6,3,7,8,12;
                      15,14,9,2,1,10,8,4,3,7,6,5,11,12,13;
                      15,14,13,8,1,12,11,7,2,10,6,3,5,4,9;
                      15,14,9,3,1,13,8,4,2,7,6,5,10,11,12;
                      15,14,10,3,1,12,9,4,2,7,6,5,8,11,13;
                      15,14,9,2,1,13,8,4,3,7,6,5,10,11,12;
                      15,13,11,9,1,12,10,8,2,7,6,3,5,4,14;
                      15,4,3,2,1,5,6,7,8,9,10,12,11,13,14;
                      15,13,4,2,1,10,6,5,3,7,8,11,9,12,14;
                      15,14,5,2,1,12,6,4,3,7,8,9,10,11,13;
                      15,14,12,9,1,13,11,8,2,10,6,3,5,4,7;
                      15,14,13,4,1,12,11,5,2,10,6,3,7,8,9;
                      15,11,10,2,1,9,8,5,3,7,6,4,12,13,14;
                      15,12,9,2,1,10,8,4,3,7,6,5,11,13,14;
                      15,12,4,2,1,8,6,5,3,7,9,11,10,13,14;
                      15,14,11,9,1,12,10,7,2,8,6,3,5,4,13;
                      15,13,9,2,1,11,8,4,3,7,6,5,10,12,14;
                      15,11,4,2,1,8,6,5,3,7,9,10,12,13,14;
                      15,14,12,3,1,13,10,4,2,8,6,5,7,9,11;
                      15,13,10,3,1,12,9,4,2,7,6,5,8,11,14;
                      15,14,13,3,1,12,11,5,2,10,6,4,7,8,9;
                      15,13,11,10,1,12,9,8,2,7,6,3,5,4,14;
                      15,14,12,4,1,13,11,5,2,9,6,3,7,8,10;
                      15,14,13,8,1,11,10,7,2,9,6,3,5,4,12;
                      15,14,11,9,1,13,10,8,2,7,6,3,5,4,12;
                      15,11,10,2,1,9,8,4,3,7,6,5,12,13,14;
                      15,12,5,2,1,9,6,4,3,7,8,10,11,13,14;
                      15,13,11,8,1,12,10,7,2,9,6,3,5,4,14;
                      15,14,11,10,1,13,9,8,2,7,6,3,5,4,12;
                      15,13,12,10,1,11,9,8,2,7,6,3,5,4,14;
                      15,14,12,4,1,13,11,5,2,10,6,3,7,8,9;
                      15,13,9,3,1,11,8,4,2,7,6,5,10,12,14;
                      15,13,5,3,1,12,6,4,2,7,8,9,10,11,14;
                      15,14,5,2,1,13,6,4,3,7,8,9,10,11,12;
                      15,14,12,3,1,13,11,4,2,9,6,5,7,8,10;
                      15,13,5,3,1,12,6,4,2,7,8,10,9,11,14;
                      15,13,11,4,1,12,10,5,2,8,6,3,7,9,14;
                      15,14,10,2,1,13,8,4,3,7,6,5,9,11,12;
                      15,13,11,4,1,12,9,5,2,8,6,3,7,10,14;
                      15,14,12,11,1,13,10,9,2,7,6,3,5,4,8;
                      15,14,12,10,1,13,11,9,2,8,6,3,5,4,7;
                      15,11,3,2,1,7,6,5,4,8,9,10,12,13,14;
                      15,14,9,3,1,12,8,4,2,7,6,5,10,11,13;
                      15,14,5,2,1,10,6,4,3,7,8,9,11,12,13;
                      15,12,3,2,1,8,6,5,4,7,9,11,10,13,14;
                      15,13,12,2,1,9,8,5,3,7,6,4,10,11,14;
                      15,13,11,9,1,12,10,7,2,8,6,3,5,4,14;
                      15,14,13,4,1,12,11,5,2,9,6,3,7,8,10;
                      15,13,11,2,1,10,8,5,3,7,6,4,9,12,14;
                      15,14,12,8,1,11,10,7,2,9,6,3,5,4,13;
                      15,14,11,8,1,13,10,7,2,9,6,3,5,4,12;
                      15,14,4,2,1,11,6,5,3,7,8,9,10,12,13;
                      15,14,13,2,1,12,10,4,3,8,6,5,7,9,11;
                      15,4,3,2,1,5,6,8,9,7,10,11,12,13,14;
                      15,14,11,4,1,13,10,5,2,8,6,3,7,9,12;
                      15,14,11,2,1,12,9,4,3,7,6,5,8,10,13;
                      15,13,12,11,1,10,9,7,2,8,6,3,5,4,14;
                      15,13,11,4,1,12,10,5,2,9,6,3,7,8,14;
                      15,14,11,3,1,13,9,4,2,7,6,5,8,10,12;
                      15,14,12,8,1,13,11,7,2,9,6,3,5,4,10;
                      15,13,3,2,1,7,6,5,4,8,9,11,10,12,14;
                      15,14,12,9,1,13,11,8,2,10,7,3,6,4,5;
                      15,13,9,2,1,10,8,4,3,7,6,5,11,12,14;
                      15,10,3,2,1,7,6,5,4,8,9,11,12,13,14;
                      15,13,4,3,1,10,6,5,2,7,8,11,9,12,14;
                      15,14,11,3,1,12,9,5,2,8,6,4,7,10,13;
                      15,14,13,2,1,11,9,5,3,8,6,4,7,10,12;
                      15,14,12,9,1,13,11,8,2,10,7,3,5,4,6;
                      15,14,4,3,1,13,6,5,2,7,8,10,9,11,12;
                      15,14,11,8,1,12,10,7,2,9,6,3,5,4,13;
                      15,14,12,3,1,11,9,5,2,8,6,4,7,10,13;
                      15,13,3,2,1,10,6,5,4,7,8,11,9,12,14;
                      15,4,3,2,1,5,6,7,8,9,10,11,12,13,14;
                      15,4,3,2,1,5,6,8,11,7,9,12,10,13,14;
                      15,11,5,2,1,9,6,4,3,7,8,10,12,13,14;
                      15,14,11,2,1,9,8,4,3,7,6,5,10,12,13;
                      15,14,13,10,1,11,9,7,2,8,6,3,5,4,12;
                      15,14,12,3,1,13,11,4,2,10,6,5,7,8,9;
                      15,14,12,10,1,13,11,9,2,7,6,3,5,4,8;
                      15,14,13,12,1,11,10,8,2,9,7,3,6,4,5;
                      15,14,4,2,1,12,6,5,3,7,8,9,10,11,13;
                      15,14,13,9,1,12,11,7,2,8,6,3,5,4,10;
                      15,14,13,11,1,12,10,7,2,9,6,3,5,4,8;
                      15,13,12,4,1,10,9,5,2,8,6,3,7,11,14;
                      15,13,12,4,1,11,9,5,2,8,6,3,7,10,14;
                      15,14,12,4,1,11,10,5,2,9,6,3,7,8,13;
                      15,13,4,2,1,9,6,5,3,7,8,10,11,12,14;
                      15,13,11,3,1,12,9,5,2,7,6,4,8,10,14;
                      15,14,12,9,1,13,11,7,2,8,6,3,5,4,10;
                      15,12,4,2,1,9,6,5,3,7,8,11,10,13,14;
                      15,14,13,9,1,12,11,8,2,10,7,3,5,4,6;
                      15,14,11,4,1,12,10,5,2,8,6,3,7,9,13;
                      15,4,3,2,1,5,6,7,9,8,10,12,11,13,14;
                      15,14,13,12,1,11,10,7,2,8,6,3,5,4,9;
                      15,11,5,2,1,8,6,4,3,7,9,10,12,13,14;
                      15,13,12,3,1,11,9,5,2,7,6,4,8,10,14;
                      15,14,13,10,1,12,11,8,2,9,6,3,5,4,7;
                      15,14,12,2,1,13,10,4,3,7,6,5,8,9,11;
                      15,14,12,11,1,10,9,8,2,7,6,3,5,4,13;
                      15,14,4,2,1,10,6,5,3,7,8,9,11,12,13;
                      15,14,11,10,1,12,9,8,2,7,6,3,5,4,13;
                      15,11,3,2,1,8,6,5,4,7,9,10,12,13,14;
                      15,13,12,4,1,11,10,5,2,8,6,3,7,9,14;
                      15,12,10,2,1,9,8,4,3,7,6,5,11,13,14;
                      15,9,3,2,1,7,6,5,4,8,10,12,11,13,14;
                      15,13,12,11,1,10,9,8,2,7,6,3,5,4,14;
                      15,14,11,3,1,12,10,5,2,8,6,4,7,9,13;
                      15,14,12,8,1,13,10,7,2,9,6,3,5,4,11;
                      15,13,12,9,1,11,10,7,2,8,6,3,5,4,14;
                      15,14,3,2,1,11,6,5,4,7,8,9,10,12,13;
                      15,14,12,2,1,13,11,4,3,10,6,5,7,8,9;
                      15,14,12,2,1,13,11,4,3,9,6,5,7,8,10;
                      15,13,12,3,1,10,9,5,2,7,6,4,8,11,14;
                      15,13,5,3,1,11,6,4,2,7,8,9,10,12,14;
                      15,14,9,2,1,12,8,4,3,7,6,5,10,11,13;
                      15,14,10,3,1,13,8,4,2,7,6,5,9,11,12;
                      15,4,3,2,1,5,6,8,10,7,9,11,12,13,14;
                      15,14,12,2,1,13,11,4,3,7,6,5,8,9,10;
                      15,14,3,2,1,7,6,5,4,8,9,10,11,12,13;
                      15,14,13,11,1,12,10,8,2,9,7,3,6,4,5;
                      15,14,13,10,1,12,11,8,2,9,7,3,6,4,5;
                      15,13,3,2,1,8,6,5,4,7,9,11,10,12,14;
                      15,14,13,12,1,10,9,7,2,8,6,3,5,4,11;
                      15,14,12,11,1,13,10,9,2,8,7,3,6,4,5;
                      15,14,4,2,1,13,6,5,3,7,8,10,9,11,12;
                      15,14,13,4,1,12,10,5,2,9,6,3,7,8,11;
                      15,14,12,10,1,13,9,8,2,7,6,3,5,4,11;
                      15,13,3,2,1,8,6,5,4,7,9,10,11,12,14;
                      15,14,13,11,1,12,9,7,2,8,6,3,5,4,10;
                      15,14,11,3,1,12,9,5,2,7,6,4,8,10,13;
                      15,13,12,2,1,10,9,5,3,8,6,4,7,11,14;
                      15,13,11,3,1,12,9,5,2,8,6,4,7,10,14;
                      15,14,11,2,1,12,8,4,3,7,6,5,9,10,13;
                      15,14,13,8,1,12,11,7,2,9,6,3,5,4,10;
                      15,14,3,2,1,9,6,5,4,7,8,10,11,12,13;
                      15,13,4,3,1,11,6,5,2,7,8,10,9,12,14;
                      15,14,13,9,1,12,11,7,2,10,6,3,5,4,8;
                      15,13,4,3,1,12,6,5,2,7,8,10,9,11,14;
                      15,14,13,2,1,12,9,4,3,7,6,5,8,10,11;
                      15,14,13,11,1,12,10,9,2,8,7,3,6,4,5;
                      15,13,10,3,1,11,8,5,2,7,6,4,9,12,14;
                      15,12,3,2,1,8,6,5,4,7,9,10,11,13,14;
                      15,14,13,10,1,12,11,8,2,9,7,3,5,4,6;
                      15,14,12,2,1,13,11,4,3,8,6,5,7,9,10;
                      15,14,13,10,1,12,11,7,2,9,6,3,5,4,8;
                      15,14,12,10,1,11,9,8,2,7,6,3,5,4,13;
                      15,14,12,10,1,13,11,8,2,7,6,3,5,4,9;
                      15,14,5,3,1,12,6,4,2,7,8,10,9,11,13;
                      15,13,3,2,1,7,6,5,4,8,9,10,11,12,14;
                      15,12,3,2,1,7,6,5,4,8,9,11,10,13,14;
                      15,14,13,2,1,12,8,4,3,7,6,5,9,10,11;
                      15,13,4,2,1,11,6,5,3,7,8,10,9,12,14;
                      15,14,12,9,1,11,10,7,2,8,6,3,5,4,13;
                      15,14,13,9,1,11,10,7,2,8,6,3,5,4,12;
                      15,14,13,12,1,11,10,9,2,7,6,3,5,4,8;
                      15,14,13,12,1,11,10,9,2,8,7,3,6,4,5;
                      15,12,3,2,1,7,6,5,4,8,9,10,11,13,14;
                      15,14,13,2,1,11,10,5,3,9,6,4,7,8,12;
                      15,14,13,2,1,12,11,4,3,10,6,5,7,8,9;
                      15,14,12,9,1,13,11,7,2,10,6,3,5,4,8;
                      15,14,12,3,1,13,11,5,2,9,6,4,7,8,10;
                      15,13,10,3,1,12,9,5,2,7,6,4,8,11,14;
                      15,14,3,2,1,7,6,5,4,8,9,11,10,12,13;
                      15,13,12,2,1,10,8,5,3,7,6,4,9,11,14;
                      15,13,10,2,1,9,8,4,3,7,6,5,11,12,14;
                      15,14,12,3,1,13,11,4,2,7,6,5,8,9,10;
                      15,14,3,2,1,12,6,5,4,7,8,10,9,11,13;
                      15,14,13,2,1,12,11,4,3,8,6,5,7,9,10;
                      15,14,13,2,1,9,8,4,3,7,6,5,10,11,12;
                      15,14,3,2,1,11,6,5,4,7,8,10,9,12,13;
                      15,10,3,2,1,7,6,5,4,8,9,12,11,13,14;
                      15,14,13,10,1,12,9,7,2,8,6,3,5,4,11;
                      15,14,11,2,1,13,9,4,3,7,6,5,8,10,12;
                      15,14,10,2,1,9,8,4,3,7,6,5,11,12,13;
                      15,14,13,2,1,11,8,4,3,7,6,5,9,10,12;
                      15,13,10,2,1,11,8,4,3,7,6,5,9,12,14;
                      15,13,11,2,1,10,8,4,3,7,6,5,9,12,14;
                      15,14,13,12,1,11,10,7,2,9,6,3,5,4,8;
                      15,14,12,3,1,13,10,5,2,9,6,4,7,8,11;
                      15,14,13,2,1,12,11,4,3,7,6,5,8,9,10;
                      15,14,12,9,1,13,10,8,2,7,6,3,5,4,11;
                      15,14,3,2,1,13,6,5,4,7,8,10,9,11,12;
                      15,14,3,2,1,8,6,5,4,7,9,10,11,12,13;
                      15,14,3,2,1,10,6,5,4,7,8,9,11,12,13;
                      15,13,11,2,1,9,8,5,3,7,6,4,10,12,14;
                      15,14,13,2,1,10,9,5,3,7,6,4,8,11,12;
                      15,14,12,3,1,13,11,5,2,10,6,4,7,8,9;
                      15,14,12,2,1,10,8,4,3,7,6,5,9,11,13;
                      15,14,4,2,1,12,6,5,3,7,8,10,9,11,13;
                      15,14,12,10,1,13,11,7,2,8,6,3,5,4,9;
                      15,14,13,9,1,12,11,8,2,10,6,3,5,4,7;
                      15,14,13,3,1,12,10,5,2,9,6,4,7,8,11;
                      15,14,12,2,1,9,8,4,3,7,6,5,10,11,13;
                      15,14,12,9,1,13,11,8,2,7,6,3,5,4,10;
                      15,14,12,10,1,13,11,8,2,9,6,3,5,4,7;
                      15,14,12,3,1,13,10,5,2,8,6,4,7,9,11;
                      15,14,4,2,1,13,6,5,3,7,8,9,10,11,12;
                      15,14,10,3,1,12,8,4,2,7,6,5,9,11,13;
                      15,13,11,2,1,9,8,4,3,7,6,5,10,12,14;
                      15,13,3,2,1,9,6,5,4,7,8,10,11,12,14;
                      15,12,4,2,1,8,6,5,3,7,9,10,11,13,14;
                      15,14,13,2,1,11,10,5,3,8,6,4,7,9,12;
                      15,13,12,3,1,10,9,5,2,8,6,4,7,11,14;
                      15,13,12,10,1,11,9,7,2,8,6,3,5,4,14;
                      15,11,3,2,1,7,6,5,4,8,9,12,10,13,14;
                      15,13,4,2,1,9,6,5,3,7,8,11,10,12,14;
                      15,14,13,8,1,12,10,7,2,9,6,3,5,4,11;
                      15,14,3,2,1,13,6,5,4,7,8,9,10,11,12;
                      15,13,4,2,1,10,6,5,3,7,8,9,11,12,14;
                      15,13,5,3,1,11,6,4,2,7,8,10,9,12,14;
                      15,14,13,3,1,11,10,5,2,9,6,4,7,8,12;
                      15,13,12,3,1,11,9,5,2,8,6,4,7,10,14;
                      15,14,12,2,1,11,9,5,3,7,6,4,8,10,13;
                      15,14,12,9,1,13,10,7,2,8,6,3,5,4,11;
                      15,14,3,2,1,9,6,5,4,7,8,11,10,12,13;
                      15,13,10,3,1,12,8,4,2,7,6,5,9,11,14;
                      15,14,13,11,1,12,10,8,2,9,6,3,5,4,7;
                      15,12,10,2,1,9,8,5,3,7,6,4,11,13,14;
                      15,14,12,2,1,9,8,5,3,7,6,4,10,11,13;
                      15,14,11,2,1,10,8,4,3,7,6,5,9,12,13;
                      15,13,10,3,1,11,8,4,2,7,6,5,9,12,14;
                      15,13,3,2,1,9,6,5,4,7,8,11,10,12,14;
                      15,14,13,2,1,10,8,4,3,7,6,5,9,11,12;
                      15,14,13,2,1,12,11,4,3,9,6,5,7,8,10;
                      15,13,12,2,1,10,9,5,3,7,6,4,8,11,14;
                      15,14,10,2,1,11,8,4,3,7,6,5,9,12,13;
                      15,14,13,2,1,12,11,5,3,9,6,4,7,8,10;
                      15,14,12,11,1,13,9,8,2,7,6,3,5,4,10;
                      15,14,13,2,1,11,9,5,3,7,6,4,8,10,12;
                      15,14,12,2,1,11,9,5,3,8,6,4,7,10,13;
                      15,14,12,10,1,11,9,7,2,8,6,3,5,4,13;
                      15,14,12,3,1,11,10,5,2,8,6,4,7,9,13;
                      15,14,12,2,1,10,8,5,3,7,6,4,9,11,13;
                      15,14,5,2,1,11,6,4,3,7,8,9,10,12,13;
                      15,14,11,9,1,13,10,7,2,8,6,3,5,4,12;
                      15,4,3,2,1,5,6,8,10,7,9,12,11,13,14;
                      15,14,10,2,1,12,8,4,3,7,6,5,9,11,13;
                      15,14,13,12,1,11,10,8,2,7,6,3,5,4,9;
                      15,14,12,4,1,11,10,5,2,8,6,3,7,9,13;
                      15,4,3,2,1,5,6,7,10,8,9,12,11,13,14;
                      15,14,3,2,1,10,6,5,4,7,8,11,9,12,13;
                      15,14,13,12,1,11,10,8,2,9,7,3,5,4,6;
                      15,14,13,12,1,10,9,8,2,7,6,3,5,4,11;
                      15,14,12,11,1,13,10,8,2,7,6,3,5,4,9;
                      15,14,12,10,1,13,11,7,2,9,6,3,5,4,8;
                      15,14,12,11,1,13,10,9,2,8,6,3,5,4,7;
                      15,14,13,2,1,11,9,4,3,7,6,5,8,10,12;
                      15,14,13,11,1,12,10,9,2,8,7,3,5,4,6;
                      15,14,13,11,1,10,9,7,2,8,6,3,5,4,12;
                      15,14,12,11,1,13,10,9,2,8,7,3,5,4,6;
                      15,14,12,10,1,13,11,8,2,9,7,3,5,4,6;
                      15,14,13,2,1,10,9,5,3,8,6,4,7,11,12;
                      15,14,12,2,1,13,9,4,3,7,6,5,8,10,11;
                      15,14,12,2,1,10,9,5,3,8,6,4,7,11,13;
                      15,13,11,3,1,10,9,5,2,7,6,4,8,12,14;
                      15,14,13,2,1,10,8,5,3,7,6,4,9,11,12;
                      15,14,13,11,1,12,10,8,2,9,7,3,5,4,6;
                      15,14,13,2,1,12,10,4,3,7,6,5,8,9,11;
                      15,14,12,2,1,11,8,4,3,7,6,5,9,10,13;
                      15,14,13,10,1,12,11,7,2,8,6,3,5,4,9;
                      15,14,13,11,1,12,9,8,2,7,6,3,5,4,10;
                      15,14,13,11,1,12,10,9,2,8,6,3,5,4,7;
                      15,14,13,11,1,10,9,8,2,7,6,3,5,4,12;
                      15,14,4,2,1,11,6,5,3,7,8,10,9,12,13;
                      15,14,13,12,1,11,9,8,2,7,6,3,5,4,10;
                      15,14,13,12,1,11,10,9,2,8,6,3,5,4,7;
                      15,14,12,2,1,10,9,5,3,7,6,4,8,11,13;
                      15,14,13,2,1,12,10,5,3,9,6,4,7,8,11;
                      15,14,12,3,1,11,9,5,2,7,6,4,8,10,13;
                      15,14,13,2,1,12,10,5,3,8,6,4,7,9,11;
                      15,14,13,10,1,12,9,8,2,7,6,3,5,4,11;
                      15,14,13,3,1,12,10,5,2,8,6,4,7,9,11;
                      15,14,13,11,1,12,10,7,2,8,6,3,5,4,9;
                      15,14,12,2,1,13,10,4,3,8,6,5,7,9,11;
                      15,14,12,11,1,10,9,7,2,8,6,3,5,4,13;
                      15,14,13,11,1,12,10,8,2,7,6,3,5,4,9;
                      15,14,13,10,1,11,9,8,2,7,6,3,5,4,12;
                      15,13,10,3,1,11,9,5,2,7,6,4,8,12,14;
                      15,14,3,2,1,8,6,5,4,7,9,11,10,12,13;
                      15,14,13,12,1,11,10,9,2,8,7,3,5,4,6;
                      15,14,13,3,1,11,10,5,2,8,6,4,7,9,12;
                      15,14,13,12,1,11,10,8,2,9,6,3,5,4,7;
                      15,14,12,2,1,11,9,4,3,7,6,5,8,10,13;
                      15,14,13,11,1,12,10,9,2,7,6,3,5,4,8;
                      15,14,3,2,1,12,6,5,4,7,8,9,10,11,13;
                      15,14,13,12,1,11,9,7,2,8,6,3,5,4,10];

% The d=1 case is trivial
elseif d == 1
    Y = cell(1);
    if mu >= 0
        Y{1} = mu;
        pos = true;
    else
        Y{1} = NaN;
        pos = false;
    end
    return

else
    err = MException('DecompEigToPSDConstruction:invalidDimension', 'The mu input to DecompEigToPSDConstruction must be a real vector of d(d+1)/2 eigenvalues, for some integer 1 <= d <= 5.');
    throw(err)
end

numPermutations = size(matricizations,1);

% p_d holds the p^\downarrow vector
p_d = zeros(1,len);
% The last entry of p^\downarrow(v) is the last entry in v
p_d(len) = mu(len);
for j=(len-1):-1:1
    % Work backwards, the previous entry is the current entry in p plus
    % the jth entry in v
    p_d(j) = p_d(j+1) + mu(j);
end

% Preallocate the cell array to hold the Y_k's
Y = cell(1,numPermutations);
for j=1:numPermutations
    Y{j} = zeros(d); % Which are d*d matrices
end

cvx_begin sdp quiet
    % X is a cvx holding chamber for the Y_k's
    variable X(d, numPermutations*d)
    % For each Y_k
    for k=1:numPermutations
        % Assign its values from the holding chamber X to the correct
        % kth place in Y
        Y{k} = X(:,(k-1)*d+1:k*d); %#ok<USENS>

        % The following performs vec = L*_k(Y_k), where the specific
        % (adjoint) matricization is chosen as the kth entry in the
        % matricization array.
        counter=1;
        % unordered_vec is a temporary vector to hold the entries
        % of Y_k until they are sorted correctly
        % For every entry ii,jj in the upper triangular part of Y_k
        for ii=1:d
            for jj=ii:d
                % Place the corresponding value into unordered_vec
                % Note: we are indexing in X, the holding chamber,
                % rather than Y so that we don't have to worry about
                % types
                unordered_vec(counter) = X(ii,(k-1)*d+jj);
                counter = counter + 1;
            end
        end
        % vec holds the output of L*_k(Y_k), where the order is
        % designated by matricizations(k,:)
        vec(matricizations(k,:)) = 2*unordered_vec;

        % p(k,:) holds the p^\uparrow vector of L*_k(Y{k})
        % The last entry of p^\uparrow(v) is v_1
        p(k,len) = vec(1);
        for j=1:len-1
            % Work backwards, the previous entry is the current entry
            % in p plus the jth entry in v
            p(k,len-j) = p(k,len-j+1) + vec(j+1);
        end

    end
subject to
    % For all k
    for k=1:numPermutations
        % Y_k must be symmetric
        Y{k} == Y{k}'
        % and Y_k must be positive semidefinite
        % Note: we double Y{k} here to ensure cvx realizes it's symmetric,
        % which avoids many warnings
        (Y{k} + Y{k}') >= 0
    end
    % p_u holds the sum of the k p^\uparrow vectors 
    if numPermutations > 1
        p_u = sum(p);
    else
        p_u = p;
    end

    for j=1:len
        % The inequalities are entrywise
        p_u(j) <= p_d(j)
    end
cvx_end

if cvx_optval == 0
    pos = true;
else
    pos = false;
    for j=1:size(Y,2)
        Y{j} = NaN;
    end
end

end

