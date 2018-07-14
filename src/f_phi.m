% phi-angle in transverse xy-plane
% 
% input: 4-momentum p = [px; py; pz; E]
% 
% mikael.mieskolainen@cern.ch, 13/07/2018

function y = phi(p)
    y = atan2(p(2), p(1));         % y, x
end
