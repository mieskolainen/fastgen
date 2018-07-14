% Calculate shannon entropy for distributions
%
% Input:    X = Histogram (1D, 2D, ...)
% Output:   S = Shannon entropy in bits
%
% mikael.mieskolainen@cern.ch, 13/07/2018

function S = shannonentropy(X)

X = X(:);
p = X / sum(X);
p = p(p > 0);
S = -sum(p.*log2(p));

end