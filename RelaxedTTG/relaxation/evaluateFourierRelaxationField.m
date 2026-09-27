function u = evaluateFourierRelaxationField(field,xk,xl)
%EVALUATEFOURIERRELAXATIONFIELD Evaluate a trigonometric relaxation field.
%
% u = Re sum_m c_m exp(i(Gk_m.xk + Gl_m.xl)) for 2-by-M inputs xk, xl.
% Evaluation is chunked so that at most field.chunkElements complex
% exponentials are held at once. gpuArray inputs are supported.
assert(size(xk,1) == 2 && isequal(size(xk),size(xl)), ...
    'xk and xl must be matching 2-by-M arrays.');
numberPoints = size(xk,2);
numberModes = size(field.Gk,2);
Gk = cast(field.Gk,'like',real(xk));
Gl = cast(field.Gl,'like',real(xk));
c = cast(field.coefficients,'like',complex(real(xk)));
u = zeros(2,numberPoints,'like',real(xk));
chunk = max(1,floor(field.chunkElements/max(numberModes,1)));
for first = 1:chunk:numberPoints
    range = first:min(first+chunk-1,numberPoints);
    phase = exp(1i*(xk(:,range).'*Gk+xl(:,range).'*Gl));
    u(:,range) = real(phase*c).';
end
end
