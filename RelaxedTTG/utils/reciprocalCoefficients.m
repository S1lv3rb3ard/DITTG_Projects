function coefficients = reciprocalCoefficients(G,B)
%RECIPROCALCOEFFICIENTS Integer n such that each row of G equals (B*n)'.
coefficientsReal = (B\G.').';
coefficients = round(coefficientsReal);
if isempty(coefficients)
    return
end
residual = max(abs(coefficientsReal-coefficients),[],'all');
assert(residual < 1e-8, ...
    'A reciprocal vector is not in the expected lattice.');
end
