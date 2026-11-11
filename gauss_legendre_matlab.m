function [omega,w,info] = gauss_legendre_matlab(n)
%GAUSS_LEGENDRE_MATLAB Base-MATLAB Gauss--Legendre rule on [0,1].
% Newton iteration in the angle variable, with a three-term recurrence.
% Only the positive roots are computed; symmetry supplies the other half.
% Reinsch's difference recurrence near x=1 avoids endpoint cancellation.
% Storage is O(n), work is O(n^2); no external library or dense eigensolve.
validateattributes(n,{'numeric'},{'scalar','integer','positive','finite'});
timer=tic;
m=floor(n/2);
theta=pi*((1:m)'-0.25)/(n+0.5);
converged=isempty(theta);
iterations=0;
for k=1:12
    if isempty(theta), break; end
    [pn,q]=legendre_value(n,theta);
    step=pn.*sin(theta)./(n*q);
    theta=theta+step;
    iterations=k;
    if all(abs(step)<=32*eps*max(theta,1/n))
        converged=true;
        break
    end
end
assert(converged,'quadrature:newton','Legendre root iteration did not converge.');
[pn,q]=legendre_value(n,theta);
lower=sin(theta/2).^2;
halfWeights=(sin(theta)./(n*q)).^2;
omega=[lower;flipud(1-lower)];
w=[halfWeights;flipud(halfWeights)];
if mod(n,2)==1
    [~,q0]=legendre_value(n,pi/2);
    omega=[lower;0.5;flipud(1-lower)];
    w=[halfWeights;1/(n*q0)^2;flipud(halfWeights)];
end
assert(numel(omega)==n && all(diff(omega)>0) && ...
    all(omega>0 & omega<1) && all(w>0 & isfinite(w)), ...
    'quadrature:invalidRule','Invalid quadrature nodes or weights.');
massError=abs(sum(w)-1);
degrees=0:40;
momentErrors=arrayfun(@(k)abs(sum(w.*omega.^k)-1/(k+1)),degrees);
assert(massError<5e-12,'quadrature:mass','Quadrature weights do not sum to one.');
% Only degrees <=2n-1 have an exactness requirement for a small rule.
maxMomentError=max(momentErrors(degrees<=2*n-1));
assert(maxMomentError<5e-12,'quadrature:moments','Polynomial moment check failed.');
info=struct('generator','gauss_legendre_matlab', ...
    'algorithm','angle Newton / Legendre recurrence / endpoint difference recurrence', ...
    'matlabVersion',version,'n',n,'iterations',iterations, ...
    'seconds',toc(timer),'massError',massError,'maxMomentError',maxMomentError, ...
    'maxPolynomialResidual',max([0;abs(pn)]));
end

function [pn,q]=legendre_value(n,theta)
% q=P_{n-1}(cos(theta))-cos(theta)*P_n(cos(theta)).
x=cos(theta);
near=theta<0.25;
ordinary=~near;
pn=zeros(size(theta)); q=pn;
if any(ordinary)
    z=x(ordinary); p0=ones(size(z)); p1=z;
    for k=2:n
        p2=((2*k-1)*z.*p1-(k-1)*p0)/k;
        p0=p1; p1=p2;
    end
    pn(ordinary)=p1;
    q(ordinary)=p0-z.*p1;
end
if any(near)
    delta=-2*sin(theta(near)/2).^2;
    p=1+delta; difference=delta;
    for k=2:n
        difference=((2*k-1)*delta.*p+(k-1)*difference)/k;
        p=p+difference;
    end
    pn(near)=p;
    q(near)=-difference-delta.*p;
end
end
