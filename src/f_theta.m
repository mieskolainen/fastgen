% theta-angle over beam-axis
% 
% input: 4-momentum p = [px; py; pz; E]
% 
% mikael.mieskolainen@cern.ch, 13/07/2018

function y = theta(p)
    y = atan2(norm(p(1:2)), p(3)); % pt, z
end
