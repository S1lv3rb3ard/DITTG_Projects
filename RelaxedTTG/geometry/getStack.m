function stack = getStack(aG,theta)
%GETSTACK Construct direct/reciprocal data for a three-layer graphene stack.
assert(isscalar(aG) && aG > 0,'aG must be positive.');
assert(isvector(theta) && numel(theta) == 3, ...
    'theta must contain the three layer angles in degrees.');

theta = theta(:).';
stack.theta = theta;
stack.aG = aG;
A_G = aG*[1,1/2;0,sqrt(3)/2];
B_G = 2*pi*inv(A_G).';
tau_B = A_G*[1/3;1/3];
K_G = B_G*[2/3;1/3];
pK_G = B_G*[1/3;2/3];

stack.V = sqrt(abs(det(B_G)));
stack.A = cell(1,3);
stack.B = cell(1,3);
stack.K = zeros(2,3);
stack.pK = zeros(2,3);
stack.tau = zeros(2,3);

for layer = 1:3
    R = rotationMatrixDegrees(theta(layer));
    stack.A{layer} = R*A_G;
    stack.B{layer} = R*B_G;
    stack.K(:,layer) = R*K_G;
    stack.pK(:,layer) = R*pK_G;
    stack.tau(:,layer) = R*tau_B;
end

scale = max(abs(stack.tau),[],'all');
stack.tau(abs(stack.tau) < 100*eps(scale)) = 0;
end
