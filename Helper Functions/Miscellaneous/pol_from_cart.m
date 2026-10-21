function POLAR = pol_from_cart(V)
% POL_FROM_CART  Magnitude, azimuth and elevation of a vector given in a
%   z-down frame such as north, east, down.
%
%   POLAR = pol_from_cart(V)
%
%   Azimuth is measured in the x-y plane from the x-axis towards the y-axis,
%   elevation from the x-y plane, positive against the z-axis (up in NED).
%   For a velocity in NED this gives speed, heading and flight path angle.
%
%   Input
%     V       [3x1]  Vector [Vx; Vy; Vz]
%
%   Output
%     POLAR   [3x1]  [magnitude; azimuth (rad); elevation (rad)]
%                    Azimuth in (-pi, pi], elevation in [-pi/2, pi/2].
%                    Azimuth is 0 for a purely vertical vector.

%   Author:    Thomas Harald Müller
%   GitHub:    https://github.com/ThomasHaraldMueller
%   Copyright: (c) 2026 Thomas Harald Müller, Hamburg, Germany
%   Free to use and modify, provided this notice is kept and the author is credited.

horiz = hypot(V(1), V(2));

POLAR = [ norm(V);
          atan2(V(2), V(1));
          atan2(-V(3), horiz) ];

end