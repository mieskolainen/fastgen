function checkrf(p, name)

% Check is a rest frame
summ = zeros(3,1);
for k = 1:length(p)
    summ = summ + p{k}(1:3);
end

if (norm(summ) > 1e-9)
    fprintf('%s: Not a rest frame %0.3f! \n', name, norm(summ));
end

end