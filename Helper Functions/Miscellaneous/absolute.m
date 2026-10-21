function ret = absolute(vec)
% ABSOLUTE  Euclidean norm of a vector.
%
%   ret = absolute(vec)
%
%   Wrapper around norm for vectors, kept for the helper struct F.
%
%   Input
%     vec   [Nx1] or [1xN]  Vector
%
%   Output
%     ret   [1x1]           Magnitude

ret = norm(vec(:));

end