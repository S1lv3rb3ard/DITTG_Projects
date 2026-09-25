function coefficient = Jackson(P)
%JACKSON Jackson damping coefficients for orders zero through P.
order = 0:P;
coefficient = ((P-order+1).*cos(pi.*order./(P+1)) ...
    +sin(pi.*order./(P+1)).*cot(pi./(P+1)))./(P+1);
coefficient(order >= P) = 0;
end
