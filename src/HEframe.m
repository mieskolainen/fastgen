% Lorentz transform to Helicity frame
%
% "Quantization axis defined by the resonance system 3-momentum vector in
% the colliding beams frame (lab frame)"
% ------------------------------------------------------------------------
% 
% input:     pf = system final state 4-vectors (cell array of (4x1) vectors)
%     direction = 1 for positive, -1 for negative beam axis orientations
%         sqrts = cms energy, e.g. 13000, in (GeV)
%
% output: pfout = transformed final states (cell array of (4x1) vectors)
%
% 4-momentum convention is p = [px,py,pz,E] = [p(1),p(2),p(3),p(4)]
%
% mikael.mieskolainen@cern.ch, 13/07/2018

function pfout = HEframe(pf, direction, sqrts)

% 1. Central system 4-momentum as a sum
X = zeros(4,1);
for k = 1:length(pf)
   X = X + pf{k}(:); 
end

% 2. Boost each final state to the system rest frame
pfout = cell(length(pf),1);
for k = 1:length(pf)
    pfout{k} = boostroutine(X, pf{k}, -1); % note minus sign
end

% 3. Boost initial state protons
mp       = 0.938; % proton mass, GeV
pb1      = [0; 0; direction * sqrt(sqrts^2/4.0 - mp^2); sqrts/2];
pb2      = [0; 0; -pb1(3); pb1(4)];
pb1boost = boostroutine(X, pb1, -1);  % note minus sign
pb2boost = boostroutine(X, pb2, -1);  % note minus sign

% 4. @@ THE POLARIZATION AXIS DEFINITION @@
zaxis = unit( -(pb1boost(1:3) + pb2boost(1:3)) );
yaxis = unit( cross( unit(pb2boost(1:3)), unit(pb1boost(1:3))) );
xaxis = unit( cross( yaxis, zaxis) ); % x = y [cross] z

% 5. Create SO(3) rotation matrix for the new coordinate axes and rotate
R = [xaxis, yaxis, zaxis]'; % Axes as columns and transpose (=inverse)
for k = 1:length(pf)
    p3new = R * pfout{k}(1:3); % Spatial part rotation Rp |-> p'
    pfout{k} = [p3new; pfout{k}(4)]; % Full 4-momentum [px; py; pz; E]
end

%}

%{
% ALTERNATIVE WAY (THE USUAL WAY, ACTUALLY)

% 2. Get Euler angles for the system
Z_angle = - f_phi(X);   % note minus sign
Y_angle = - f_theta(X); % note minus sign

% 3. Rotate final states
pfout = pf;
for k = 1:length(pf)
    pfout{k} = rotateXYZ(pfout{k}, Z_angle, 3);  % z-rotation
    pfout{k} = rotateXYZ(pfout{k}, Y_angle, 2);  % y-rotation
end

% 4. New central system 4-momentum in the rotated frame
XNEW = zeros(4,1);
for k = 1:length(pfout)
   XNEW = XNEW + pfout{k}(:); 
end

% This should give [0 0 +pz E], that is, we have rotated the frame so that
% the system is aligned on the z-axis (CHECK for DEBUG)
% XNEW

% 5. Boost each final state to the central system rest frame
for k = 1:length(pfout)
    pfout{k} = boostroutine(XNEW, pfout{k}, -1); % note minus sign
end
%}

checkrf(pfout, 'HEframe'); % Check numerically that we have a rest frame
end

% Return unit vector
function y = unit(x)
    y = x / norm(x);
end