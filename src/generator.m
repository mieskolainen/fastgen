% Fast/toy MC generator for exclusive systems (no proton legs) + 2-body decay
% ------------------------------------------------------------------------
%
% No flat lorentz invariant measure dLIPS by construction enforced here,
% that is left for the user to keep track of / modify. See PDG/kinematics.
% Kinematics conserves 4-momentum.
 
% Note duality exponential in X^2 <=> Rayleigh distributed in X.
%
% 
% Input:    ptmode = 1,2,3 (flat in pt^2, exponential in pt^2, flat in pt)
%         massmode = 1,2,3 (flat in m^2,  exponential in m^2,  flat in m)
%           lambda = system pt distribution parameter (1/<pt^2>) (mode = 2)
%            kappa = mass distribution parameter (1/<m^2>)       (mode = 2)
%           limits = sampling limits (struct)
%                    .ptmin for system transverse momentum
%                    .ptmax
%                    .mmin  for system mass
%                    .mmax
%                    .ymin  for system rapidity
%                    .ymax
%             mdec = decay daughter masses (2 x 1) vector
%
% Output:        p = System 4-momentum
%               p1 = Daughter 1 4-momentum
%               p2 = Daughter 2 4-momentum
%
% 4-momentum convention is p = [px,py,pz,E] = [p(1),p(2),p(3),p(4)]
%
% mikael.mieskolainen@cern.ch, 13/07/2018

function [p, p1, p2] = generator(ptmode, massmode, lambda, kappa, limits, mdec)

% Inline function for uniform random numbers from [a,b]
U = @(a,b) a + (b-a)*rand(1);

% Minkowski metric
g = [-1  0  0  0;
     0  -1  0  0;
     0   0 -1  0
     0   0  0  1];

% Pick Random variables --------------------------------------------------


% MASS

% Flat in m^2
if (massmode == 1)
    m2 = U(limits.mmin^2, limits.mmax^2);
    
% Exponential in m^2
elseif (massmode == 2)
    % https://en.wikipedia.org/wiki/Inverse_transform_sampling
    while (true)
        m2 = -1/kappa * log(1-rand(1));
        
        if (sqrt(m2) > limits.mmin && sqrt(m2) < limits.mmax)
            break;
        end
    end
    
% Flat in m
elseif (massmode == 3)
    m2 = U(limits.mmin, limits.mmax)^2;
end


% PT

% Flat in pt^2
if (ptmode == 1)
    pt2 = U(limits.ptmin^2, limits.ptmax^2);
    
% Exponential in pt^2
elseif (ptmode == 2)
    % https://en.wikipedia.org/wiki/Inverse_transform_sampling
    while (true)
        pt2 = -1/lambda * log(1-rand(1));
        
        if (sqrt(pt2) > limits.ptmin && sqrt(pt2) < limits.ptmax)
            break;
        end
    end
    
% Flat in pt
elseif (ptmode == 3)
    pt2 = U(limits.ptmin, limits.ptmax)^2;
end

% Pick the system rapidity flat
y = U(limits.ymin, limits.ymax);

% Pick random phi from [0,2pi]
phi = U(0, 2*pi);


% System Kinematics ------------------------------------------------------

% Transverse mass -> Hyperbolic space "hypotenuse", see e.g. PDG/kinematics
mT = sqrt(m2 + pt2);

E  = mT * cosh(y); % Energy
pz = mT * sinh(y); % Longitudinal momentum

% Transverse momentum components -> Euclidean xy-space
px = sqrt(pt2) * cos(phi);
py = sqrt(pt2) * sin(phi);

% System 4-momentum finally
p = [px; py; pz; E];


% Decay Kinematics -------------------------------------------------------

[p1,p2] = twobody(p, mdec(1), mdec(2));


% ========================================================================
% Check 4-momentum conservation

% Check invariant mass squared
m2calc = abs( m2 - p'*g*p );

% Check 4-momentum is conserved in the decay
pcalc  = sum(abs( p - (p1 + p2) ));

if (m2calc > 1e-9 || pcalc > 1e-9)
    fprintf('Warning: Event kinematics not ok! m2 = %0.2E, pnorm = %0.2E \n', ...
        m2 - m2calc, pcalc);
end
% ========================================================================

end
