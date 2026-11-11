function sol=mdm_nodes(red,times,h)
%MDM_NODES Modified Davison-Maki method for the projected NDRE.
% Restarted Radon step [U;V]=expm(h*H)*[I;Y], Y<-V/U with
% H=[D,-S;Q,-A] (Davison and Maki 1973; Angelova et al. 2020, EBA-exp;
% Behr 2021, Algorithm 6.3). Each output interval is split into equal
% substeps no longer than h; step exponentials are cached by length.
[k,ell]=size(red.Q); H=[red.D,-red.S;red.Q,-red.A];
if isfield(red,'Y0'), Y=red.Y0; else, Y=zeros(k,ell); end
history=zeros(k,ell,numel(times)); history(:,:,1)=Y;
keys=zeros(0,1); cache={}; steps=0;
for j=2:numel(times)
    dt=times(j)-times(j-1); m=max(1,ceil(dt/h-1e-12)); step=dt/m;
    hit=find(abs(keys-step)<=1e-13*step,1);
    if isempty(hit), keys(end+1,1)=step; cache{end+1}=expm(step*H); hit=numel(keys); end %#ok<AGROW>
    E=cache{hit}; E11=E(1:ell,1:ell); E12=E(1:ell,ell+1:end);
    E21=E(ell+1:end,1:ell); E22=E(ell+1:end,ell+1:end);
    for i=1:m
        U=E11+E12*Y; Y=(E21+E22*Y)/U;
    end
    steps=steps+m; history(:,:,j)=Y;
end
sol=struct('Yhistory',history,'Y',Y,'time',times,'steps',steps, ...
    'distinctExponentials',numel(keys));
end
