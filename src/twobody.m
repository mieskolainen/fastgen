% 2-body isotropic decay routine according to dLIPS_2 flat phase space
% ------------------------------------------------------------------------
%
% Input:       pmot = 4-momentum of the mother (in the lab frame, usually)
%                m1 = Daughter 1 invariant mass
%                m2 = Daughter 2 invariant mass
%
% Output:     p1,p2 = Daughter 4-momentum
%                     in the the same frame as the mother was defined
%
% 4-momentum convention is p = [px,py,pz,E] = [p(1),p(2),p(3),p(4)]
%
% mikael.mieskolainen@cern.ch, 13/07/2018

function [p1, p2] = twobody(pmot, m1, m2)

% Mother mass
m0 = sqrt(pmot(4)^2 - norm(pmot(1:3))^2);

% Energies, and momentum absolute (back to back)
e1   = 0.5 * (m0^2 + m1^2 - m2^2) / m0;
e2   = 0.5 * (m0^2 + m2^2 - m1^2) / m0;
pabs = 0.5 * sqrt( (m0 - m1 - m2) * (m0 + m1 + m2) ...
           * (m0 + m1 - m2) * (m0 - m1 + m2) ) / m0;

% Isotropic angles in a spherical system
costheta = 2.0 * rand(1) - 1.0; % [-1,1]
sintheta = sqrt(1.0 - costheta^2);
phi      = 2.0 * pi * rand(1);  % [0,2pi]

% To cartesian
pX       = pabs * sintheta * cos(phi);
pY       = pabs * sintheta * sin(phi);
pZ       = pabs * costheta;

% 4-momenta now defined in the mother rest frame
p1 =  [ pX,  pY,  pZ, e1]';
p2 =  [-pX, -pY, -pZ, e2]';

% Then boost daughters into the original frame
sign = 1;
p1 = boostroutine(pmot, p1, sign);
p2 = boostroutine(pmot, p2, sign);

end