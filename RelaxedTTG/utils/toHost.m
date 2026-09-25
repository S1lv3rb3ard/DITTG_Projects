function y = toHost(x,useGPU)
%TOHOST Gather an array only when GPU execution is enabled.
if useGPU
    y = gather(x);
else
    y = x;
end
end
