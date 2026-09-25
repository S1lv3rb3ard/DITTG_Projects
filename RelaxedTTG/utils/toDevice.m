function y = toDevice(x,useGPU)
%TODEVICE Move an array to the GPU only when requested.
if useGPU
    y = gpuArray(x);
else
    y = x;
end
end
