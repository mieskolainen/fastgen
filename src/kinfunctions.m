% mikael.mieskolainen@cern.ch, 110718

% Kinematic inline functions
f_eta = @(p) atanh( p(3) / norm(p(1:3)) );          % atanh(pz/|\vec{p}|)
f_rap = @(p) 0.5 * log((p(4)+p(3)) / (p(4)-p(3)) ); % 0.5 ln((E+pz)/(E-pz))
f_m   = @(p) sqrt(p(4)^2 - norm(p(1:3))^2);         % E^2 = m^2 + |p|^2
f_pt  = @(p) norm(p(1:2));                          % (px^2 + py^2)^{-1}
