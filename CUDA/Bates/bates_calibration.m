function sol = bates_calibration(stockprice,strike,tau,typed,price)
	%
	%
	n = length(stockprice);
	pp = gpuArray(ones(1,n));
	kernel = parallel.gpu.CUDAKernel('bates.ptx','bates.cu','kernelbates');
	kernel.GridSize = [1024 1 1];
	kernel.ThreadBlockSize = [512 1 1];
	%
	%
	function z = objf(x)
		kappa  = x(1);
		theta  = x(2);
		sigma  = x(3);
		rho    = x(4);
		v0     = x(5);
		lambda = x(6);
	        muj    = x(7);
		sigmaj = x(8);	
		%
		pp = feval(kernel,pp,stockprice,strike,tau,typed,n,0.02,kappa,theta,sigma,rho,v0,lambda,muj,sigmaj);
		pp1 = gather(pp);
		z = sum(abs(pp1 - price)./price);
	end
	%
	%
	function [c,ceq] = nlc(x)
             c = x(3)*x(3) - 2*x(1)*x(2);
             ceq = [];
         end
    	%
    	%
    	lb = [  0.1  0.001    0.001   -1.0   0.001   0.001   -1.0  0.001];
    	ub = [100.0  3.0      3.0      0.0   3.0    30.0      1.0  2.000];
    	%
    	options = optimoptions('fmincon','Display','iter', ...
                            'MaxFunctionEvaluations',10000, ...
                            'Algorithm','interior-point');
    	%
    	x0 = [20.0    0.05   0.05    -0.5   0.02  1.0 -0.1  0.02];
    	sol  = fmincon(@objf,x0,[],[],[],[],lb,ub,@nlc,options);
end
