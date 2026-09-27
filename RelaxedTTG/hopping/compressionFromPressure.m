function epsilon = compressionFromPressure(pressure)
%COMPRESSIONFROMPRESSURE Relative interlayer distance under uniaxial pressure.
%
% epsilon = compressionFromPressure(P) inverts the fit
%   P = A (exp(-B epsilon) - 1),  A = 5.73 GPa, B = 9.54,
% of Carr, Fang, Jarillo-Herrero & Kaxiras, PRB 98, 085144 (2018), for
% bilayer graphene encapsulated in hBN. P is in GPa; epsilon = d/d0 - 1 is
% negative under compression (P = 9.2 GPa gives epsilon ~ -0.10).
A = 5.73;
B = 9.54;
assert(all(pressure(:) > -A),'The pressure must exceed -%.2f GPa.',A);
epsilon = -log(1+pressure/A)/B;
end
