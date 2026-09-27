function julia = readJuliaRelaxationFiles(minimizerFile,dataFile)
%READJULIARELAXATIONFILES Read the raw binary output of example2.jl.
%
% example2.jl calls Base.write, which stores raw little-endian bytes:
%   data file:      N (Int64), theta1, theta3, E (2x2), P (2), K, G (Float64)
%   minimizer file: u (Float64, 2 x N x N x N x N x 3, column major)
% Julia and MATLAB are both column major, so reshape reproduces the layout
% u(component, s, t, v, w, layer).
fid = fopen(dataFile,'r','ieee-le');
assert(fid > 0,'Cannot open %s.',dataFile);
cleanup = onCleanup(@() fclose(fid));
julia.N = double(fread(fid,1,'int64'));
angles = fread(fid,2,'double');
julia.theta1 = angles(1);
julia.theta3 = angles(2);
julia.E = reshape(fread(fid,4,'double'),2,2);
julia.P = fread(fid,2,'double');
moduli = fread(fid,2,'double');
if numel(moduli) == 2
    julia.K = moduli(1);
    julia.G = moduli(2);
else
    julia.K = NaN;
    julia.G = NaN;
end
clear cleanup

N = julia.N;
assert(N >= 1 && N == fix(N) && N < 1e4, ...
    'Unexpected hull size N = %g in %s.',N,dataFile);
fid = fopen(minimizerFile,'r','ieee-le');
assert(fid > 0,'Cannot open %s.',minimizerFile);
cleanup = onCleanup(@() fclose(fid));
u = fread(fid,inf,'double');
clear cleanup
assert(numel(u) == 2*N^4*3, ...
    ['%s holds %d values; expected 2*N^4*3 = %d for N = %d. Check that ' ...
     'the data and minimizer files belong to the same run.'], ...
    minimizerFile,numel(u),2*N^4*3,N);
julia.u = reshape(u,[2,N,N,N,N,3]);
end
