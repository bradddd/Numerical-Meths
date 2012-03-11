%use SOR to solve 2-d backward euler diffusion equation
% x(k+1) = (I-inv(D))*x(k) + inv(D)*b
%
%
hold off;
axis auto;
C=4;
w=1.75;
n=200;
%z=linspace(-2,2,n);
%[xx yy] = meshgrid(z,z);
MAX_ITER=1000;
nsteps = 1000;
%xold = exp(-2*xx.^2).*exp(-2*yy.^2);
xold = zeros(n,n);
xold(n/2,n/2)=2.0;
forcing = xold;
xnew = zeros(n,n);

for step=1:nsteps
    step
    x = xold;
    niter=0;
    for m=1:MAX_ITER
        lastx = x;
        for i=2:n-1
            for j=2:n-1
                x(i,j) = (1-w)*x(i,j) + w*C/(4*C+1)*( x(i-1,j) + x(i+1,j) + ...
                    + x(i,j-1) + x(i,j+1)) + w/(4*C+1)*(xold(i,j)+forcing(i,j));
            end
        end
        if (mean(abs(x-lastx)) < 1.e-6)
            break;
        end
    end
    %sprintf('iterations: %d, \n', m)
    mesh(x);
%    surfc(x);
    axis manual;
    axis([0 n 0 n 0 1]);
    drawnow;
    %if (n==1)
    %    hold on;
    %end
    
    xold=x;
end
mesh(x);