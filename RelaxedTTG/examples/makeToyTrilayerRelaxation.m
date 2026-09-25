function relaxationFields = makeToyTrilayerRelaxation(stack,amplitude)
%MAKETOYTRILAYERRELAXATION Smooth periodic fields for testing the pipeline.
if nargin < 2 || isempty(amplitude)
    amplitude = 0.01*stack.aG;
end

relaxationFields = cell(1,3);
for layer = 1:3
    others = setdiff(1:3,layer,'stable');
    Ak = stack.A{others(1)};
    Al = stack.A{others(2)};
    layerPhase = 2*pi*(layer-1)/3;
    relaxationFields{layer} = @(xk,xl) ...
        toyField(xk,xl,Ak,Al,amplitude,layerPhase);
end
end

function u = toyField(xk,xl,Ak,Al,amplitude,layerPhase)
sk = mod(Ak\xk,1);
sl = mod(Al\xl,1);
u = zeros(2,size(xk,2),'like',xk);
u(1,:) = amplitude*( ...
    sin(2*pi*sk(1,:)+layerPhase) ...
    +0.60*sin(2*pi*sl(1,:)-layerPhase) ...
    +0.25*sin(2*pi*(sk(1,:)+sl(2,:))));
u(2,:) = amplitude*( ...
    sin(2*pi*sk(2,:)+layerPhase) ...
    -0.60*sin(2*pi*sl(2,:)-layerPhase) ...
    +0.25*sin(2*pi*(sk(2,:)-sl(1,:))));
end
