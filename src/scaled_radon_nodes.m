function sol=scaled_radon_nodes(red,times,scheme)
%SCALED_RADON_NODES Scaled quotient maps with optional autonomous composition.
% Both direct and composed implementations preserve identity increments.
if nargin<3, scheme="composed"; end
assert(any(scheme==["direct","composed"]),'radon:scheme','Unknown scheme.');
assert(all(isfinite(times)) && all(times>=0) && all(diff(times)>0));
[k,ell]=size(red.Q); Y0=zeros(k,ell);
if isfield(red,'Y0'), Y0=red.Y0; end
history=zeros(k,ell,numel(times)); Y=Y0; previousTime=0;
cache=containers.Map('KeyType','char','ValueType','any'); maxPerturbation=0;
for j=1:numel(times)
    if times(j)==0, history(:,:,j)=Y0; continue; end
    if scheme=="direct", delta=times(j); Y=Y0;
    else, delta=times(j)-previousTime; end
    key=sprintf('%.13e',delta);
    if isKey(cache,key), map=cache(key);
    else, map=step_map(red,delta); map.delta=delta; cache(key)=map; end
    maxPerturbation=max(maxPerturbation,abs(delta-map.delta));
    M=eye(ell)-map.S*Y;
    assert(rcond(M)>1e-14,'radon:denominator','Quotient denominator failed.');
    Y=map.Q+(map.L*Y/M)*map.R; history(:,:,j)=Y;
    previousTime=times(j);
end
sol=struct('Yhistory',history,'Y',history(:,:,end),'time',times, ...
    'scheme',scheme,'maximumCachedStepPerturbation',maxPerturbation);
end

function map=step_map(red,h)
k=size(red.A,1); ell=size(red.D,1); H=[red.D,-red.S;red.Q,-red.A];
s=max(0,ceil(log2(max(h*norm(H,1),realmin)/.25)));
B=(h/2^s)*H; term=B; delta=B;
for q=2:40
    term=term*B/q; delta=delta+term;
    if norm(term,1)<=eps*max(norm(delta,1),realmin), break; end
end
U=eye(ell)+delta(1:ell,1:ell);
dR=-(U\delta(1:ell,1:ell)); S=-(U\delta(1:ell,ell+1:end));
Q=delta(ell+1:end,1:ell)/U;
dL=delta(ell+1:end,ell+1:end)-Q*delta(1:ell,ell+1:end);
for q=1:s
    L=eye(k)+dL; R=eye(ell)+dR; U=eye(k)-Q*S; V=eye(ell)-S*Q;
    assert(min(rcond(U),rcond(V))>1e-14,'radon:doubling','Quotient doubling failed.');
    nextQ=Q+L*(U\Q)*R; nextS=S+R*(V\S)*L;
    dL=2*dL+dL*dL+L*(U\(Q*S))*L;
    dR=2*dR+dR*dR+R*(V\(S*Q))*R;
    Q=nextQ; S=nextS;
end
map=struct('Q',Q,'S',S,'L',eye(k)+dL,'R',eye(ell)+dR);
end
