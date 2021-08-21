%%  YEigIneqD3
%   Determines whether or not there exists a positive semidefinite matrix Y
%   such that several inequalities are obeyed, as specified in Theorem 5 of
%   the paper.
%   
%   This function has one required argument:
%     mu: a row or column vector of 6 real eigenvalues
%
%   pos = YEigIneqD3(mu) outputs true if there exists such a matrix Y, and
%   false otherwise.
%
%   [pos, Y] = YEigIneqD3(mu) makes the same calculation but additionally
%   outputs the matrix Y if it exists.

%   requires: cvx (http://cvxr.com/cvx/)
%   author: Logan Pipes (ldpipes@mta.ca)
%   last updated: August 20, 2021

%#ok<*NOPRT>

function [pos, Y] = YEigIneqD3(mu)

mu = sort(mu(:)', 'descend'); %#ok<UDIM>

if ~all(size(mu) == [1 6])
    err = MException('YEigIneqD3:invalidDimension', 'The mu input to YEigIneqD3 must be a vector of 6 real eigenvalues.');
    throw(err)
end

cvx_begin sdp quiet
    variable Y(3,3) symmetric semidefinite
subject to
    Y(1,1) + Y(3,3) + Y(2,3) + Y(2,2) + Y(1,2) + Y(1,3) <= mu(6) + mu(5) + mu(4) + mu(3) + mu(2) + mu(1) %#ok<NODEF>
             Y(3,3) + Y(2,3) + Y(2,2) + Y(1,2) + Y(1,3) <= mu(6) + mu(5) + mu(4) + mu(3) + mu(2)
                      Y(2,3) + Y(2,2) + Y(1,2) + Y(1,3) <= mu(6) + mu(5) + mu(4) + mu(3)
                               Y(2,2) + Y(1,2) + Y(1,3) <= mu(6) + mu(5) + mu(4)
                                        Y(1,2) + Y(1,3) <= mu(6) + mu(5)
                                                 Y(1,3) <= mu(6)
cvx_end

if cvx_optval == 0
    pos = true;
else
    pos = false;
    Y = NaN;
end

end

