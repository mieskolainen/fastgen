% Lorentz transform to "Anti-Helicity" frame
%
% "Quantization z-axis defined by the bisector vector between initial state
% p1 and (p2) (POSITIVE) directions in the (resonance) system rest frame,
% where p1 and p2 are the initial state proton 3-momentum"
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

function pfout = BEframe(pf, direction, sqrts)

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
zaxis = unit( unit(pb1boost(1:3)) + unit(pb2boost(1:3)) );
yaxis = unit( cross( unit(pb2boost(1:3)), unit(pb1boost(1:3))) );
xaxis = unit( cross( yaxis, zaxis) ); % x = y [cross] z

% 5. Create SO(3) rotation matrix for the new coordinate axes and rotate
R = [xaxis, yaxis, zaxis]'; % Axes as columns and transpose (=inverse)
for k = 1:length(pf)
    p3new = R * pfout{k}(1:3); % Spatial part rotation Rp |-> p'
    pfout{k} = [p3new; pfout{k}(4)]; % Full 4-momentum [px; py; pz; E]
end

checkrf(pfout, 'BEframe'); % Check numerically that we have a rest frame
end

% Return unit vector
function y = unit(x)
    y = x / norm(x);
end

% Bisector vector
function b = bisector(u,v)
    b = norm(u)*v + norm(v)*u;
end