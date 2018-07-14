% Lorentz Transform to Gottfried-Jackson (pseudo) frame
%
% "Quantization axis defined by the initial state proton
%  p1 (or p2) 3-momentum vector in the resonance system rest frame"
%
% "Pseudo" means here that we distinquish for historical reasons
% the "original Gottfried-Jackson" / "true Gottfried-Jackson" frame,
% where we replace proton 1 (or 2) by the true initial states such as
% fusing quarks, gluons, Pomerons etc. Clearly, access to that frame is
% only possible in fully exclusive processes / measurements.
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

function pfout = GJframe(pf, direction, sqrts)

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
zaxis = unit( pb1boost(1:3) );
yaxis = unit( cross( unit(pb2boost(1:3)), unit(pb1boost(1:3))) );
xaxis = unit( cross( yaxis, zaxis) ); % x = y [cross] z

% 5. Create SO(3) rotation matrix for the new coordinate axes and rotate
R = [xaxis, yaxis, zaxis]'; % Axes as columns and transpose (=inverse)
for k = 1:length(pf)
    p3new = R * pfout{k}(1:3); % Spatial part rotation Rp |-> p'
    pfout{k} = [p3new; pfout{k}(4)]; % Full 4-momentum [px; py; pz; E]
end

%{
% ALTERNATIVE WAY

% 4. Get rotation angles
Z_angle = - f_phi(pb1boost);
Y_angle = - f_theta(pb1boost);

% Test rotation, should give [0; 0; +pz; E] (only for debug)
% pb1boost = rotateXYZ(pb1boost, Z_angle, 3);
% pb1boost = rotateXYZ(pb1boost, Y_angle, 2)

% 5. Rotate final states
for k = 1:length(pf)
    pfout{k} = rotateXYZ(pfout{k}, Z_angle, 3);  % z-rotation
    pfout{k} = rotateXYZ(pfout{k}, Y_angle, 2);  % y-rotation
end
%}

checkrf(pfout, 'GJframe'); % Check numerically that we have a rest frame
end

% Return unit vector
function y = unit(x)
    y = x / norm(x);
end