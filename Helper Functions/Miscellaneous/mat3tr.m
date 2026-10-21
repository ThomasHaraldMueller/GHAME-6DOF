function TBD = mat3tr(psi, tht, phi)
% MAT3TR  Transformation matrix from Euler angles in the 3-2-1 sequence.
%
%   TBD = mat3tr(psi, tht, phi)
%
%   Rotates first by yaw about z, then by pitch about the new y, then by roll
%   about the new x. Applied to the geodetic frame this gives the matrix from
%   geodetic to body axes, v_B = TBD * v_D. The function euler is its inverse.
%
%   Inputs
%     psi   [1x1]  Yaw angle (rad)
%     tht   [1x1]  Pitch angle (rad)
%     phi   [1x1]  Roll angle (rad)
%
%   Output
%     TBD   [3x3]  Transformation matrix, geodetic to body

%   Author:    Thomas Harald Müller
%   GitHub:    https://github.com/ThomasHaraldMueller
%   Copyright: (c) 2026 Thomas Harald Müller, Hamburg, Germany
%   Free to use and modify, provided this notice is kept and the author is credited.

spsi = sin(psi);   cpsi = cos(psi);
stht = sin(tht);   ctht = cos(tht);
sphi = sin(phi);   cphi = cos(phi);

TBD = [ cpsi*ctht,                    spsi*ctht,                    -stht;
        cpsi*stht*sphi - spsi*cphi,   spsi*stht*sphi + cpsi*cphi,    ctht*sphi;
        cpsi*stht*cphi + spsi*sphi,   spsi*stht*cphi - cpsi*sphi,    ctht*cphi ];

end