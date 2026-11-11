function sol=incremental_centered_nodes(red,prep,times,scheme)
%INCREMENTAL_CENTERED_NODES Classical centered flow evaluated via increments.
% Kernel: a Taylor series on a small interval followed by K(2h)=K+R*K*P.
% For arbitrary initial values, expand the centered quotient about Y0.
% No subtraction L-R*L*P or Ystar+E is used. This is not a new exact flow.
if nargin<4, scheme="direct"; end
assert(any(scheme==["direct","composed"]),'incremental:scheme','Unknown scheme.');
assert(all(isfinite(times)) && all(times>=0) && all(diff(times)>0), ...
    'incremental:times','Times must be finite, nonnegative and increasing.');
[k,ell]=size(red.Q); history=zeros(k,ell,numel(times));
Y0=zeros(k,ell); if isfield(red,'Y0'), Y0=red.Y0; end
E0=Y0-prep.Ystar; Ac=prep.Ac; Dc=prep.Dc;
dP=zeros(k); dR=zeros(ell); K=zeros(ell,k); previousTime=0;
cache=containers.Map('KeyType','char','ValueType','any'); maximumTimePerturbation=0;
for j=1:numel(times)
    t=times(j);
    if t==0, history(:,:,j)=Y0; continue; end
    if scheme=="direct"
        [dP,dR,K]=kernel_step(Ac,Dc,red.S,t);
    else
        delta=t-previousTime; key=sprintf('%.13e',delta);
        if isKey(cache,key), step=cache(key);
        else
            [step.dP,step.dR,step.K]=kernel_step(Ac,Dc,red.S,delta);
            step.delta=delta; cache(key)=step;
        end
        maximumTimePerturbation=max(maximumTimePerturbation,abs(delta-step.delta));
        P=eye(k)+dP; R=eye(ell)+dR;
        K=K+R*step.K*P;
        dP=dP+step.dP+dP*step.dP; dR=dR+step.dR+dR*step.dR;
        previousTime=t;
    end
    P=eye(k)+dP; R=eye(ell)+dR; M=eye(ell)-K*E0;
    assert(rcond(M)>1e-14,'pilot:denominator','Incremental denominator failed.');
    history(:,:,j)=Y0+dP*E0+E0*dR+dP*E0*dR+P*E0*(M\(K*E0))*R;
end
sol=struct('Yhistory',history,'Y',history(:,:,end),'time',times, ...
    'scheme',scheme,'maximumCachedStepPerturbation',maximumTimePerturbation);
end

function [dP,dR,K]=kernel_step(Ac,Dc,S,t)
k=size(Ac,1); ell=size(Dc,1); scale=norm(Ac,1)+norm(Dc,1);
squarings=max(0,ceil(log2(max(t*scale,realmin)/.25))); h=t/2^squarings;
dP=exponential_increment(-h*Ac); dR=exponential_increment(-h*Dc);
term=h*S; K=term;
for q=1:40
    term=-h*(Dc*term+term*Ac)/(q+1); K=K+term;
    if norm(term,1)<=eps*max(norm(K,1),realmin), break; end
end
for q=1:squarings
    P=eye(k)+dP; R=eye(ell)+dR; K=K+R*K*P;
    dP=2*dP+dP*dP; dR=2*dR+dR*dR;
end
end

function dE=exponential_increment(B)
term=B; dE=B;
for q=2:40
    term=term*B/q; dE=dE+term;
    if norm(term,1)<=eps*max(norm(dE,1),realmin), break; end
end
end
