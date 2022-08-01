%
stockprice = gpuArray(load("stockprice.txt")');
strike = gpuArray(load("strike.txt")');
tau = gpuArray(load("tau.txt")');
typed = gpuArray(int32(load("type.txt"))');
price = load("price.txt")';
n = length(price);
%
% set kernel
% 
%kernel = parallel.gpu.CUDAKernel('heston.ptx','heston.cu','kernelheston');
%kernel.GridSize = [1024 1 1];
%kernel.ThreadBlockSize = [512 1 1];
%
%pp = gpuArray(ones(1,n));
