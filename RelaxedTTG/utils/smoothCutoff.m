function value = smoothCutoff(radius,innerRadius,outerRadius)
%SMOOTHCUTOFF Smoothly decrease from one to zero on [inner,outer].
assert(innerRadius >= 0 && outerRadius > innerRadius, ...
    'The cutoff radii must satisfy 0 <= innerRadius < outerRadius.');

value = ones(size(radius),'like',radius);
value(radius >= outerRadius) = 0;
transition = radius > innerRadius & radius < outerRadius;
t = (radius(transition)-innerRadius)/(outerRadius-innerRadius);
left = exp(-1./max(1-t,eps));
right = exp(-1./max(t,eps));
value(transition) = left./(left+right);
end
