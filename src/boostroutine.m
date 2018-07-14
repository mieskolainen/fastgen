% Lorentz boost
% ------------------------------------------------------------------------
%
% Input:   boostvec = Boost 4-momentum (4 x 1), e.g., mother 4-momentum
%                 p = 4-vector / 4-momentum to be boosted (4 x 1)
%              sign = 1 or -1 (direction of the boost, in or out)
% 
% Output       pout = Boosted vector (4 x 1)
%
% 4-momentum convention is p = [px,py,pz,E] = [p(1),p(2),p(3),p(4)]
%
% mikael.mieskolainen@cern.ch, 13/07/2018

function pout = boostroutine(boostvec, p, sign)

% Mother mass
mass = sqrt(boostvec(4)^2 - norm(boostvec(1:3))^2);

% Beta and gamma factors    
betaX = sign * boostvec(1) / boostvec(4); % px / E
betaY = sign * boostvec(2) / boostvec(4); % py / E
betaZ = sign * boostvec(3) / boostvec(4); % pz / E
gamma = boostvec(4) / mass;               % E / m

% Momentum and energy product
aux1 = betaX * p(1) + betaY * p(2) + betaZ * p(3);
aux2 = gamma * (gamma * aux1 / (1.0 + gamma) + p(4));

% Lorentz boost
pout    = zeros(4,1);
pout(1) = p(1) + aux2 * betaX;
pout(2) = p(2) + aux2 * betaY;
pout(3) = p(3) + aux2 * betaZ;
pout(4) = gamma * (p(4) + aux1);

end