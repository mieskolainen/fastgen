% Active SO(3) rotation around X, Y or Z axis
% ------------------------------------------------------------------------
%
% Input:         p  =  4-vector to be rotated
%            angle  =  Rotation angle (in radians)
%              axis =  1,2,3 (for X,Y,Z)
%
% Output:        k  =  Rotated 4-vector
% 
% 4-momentum convention is p = [px,py,pz,E] = [p(1),p(2),p(3),p(4)]
%
% mikael.mieskolainen@cern.ch, 13/07/2018

function k = rotateXYZ(p, angle, axis)

p = p(:);
c = cos(angle);
s = sin(angle);

% Rotations
if     (axis == 1) % X
    R = [1 0  0;
         0 c -s;
         0 s  c];

elseif (axis == 2) % Y
    R = [c 0 s;
         0 1 0;
        -s 0 c];

elseif (axis == 3) % Z
    R = [c -s 0;
         s  c 0;
         0  0 1];
end

% Map
k      = zeros(4,1);
k(1:3) = R * p(1:3);  % Spatial part
k(4)   = p(4);        % Energy

end