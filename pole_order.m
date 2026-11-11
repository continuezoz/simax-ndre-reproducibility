function s=pole_order(s,order)
%POLE_ORDER Reorder a fixed pole set so that every prefix samples its range.
% "monotone" keeps the given order, "leja" uses discrete Leja points of
% log|s| starting at s(1), and "interleave" uses bit-reversed indices.
n=numel(s);
switch string(order)
    case "monotone"
    case "leja"
        x=log(abs(s)); idx=zeros(1,n); idx(1)=1; used=false(1,n); used(1)=true;
        logdist=log(abs(x-x(1))+realmin);
        for j=2:n
            logdist(used)=-Inf; [~,k]=max(logdist); idx(j)=k; used(k)=true;
            logdist=logdist+log(abs(x-x(k))+realmin);
        end
        s=s(idx);
    case "interleave"
        bits=ceil(log2(n)); r=zeros(1,2^bits);
        for j=0:2^bits-1, r(j+1)=bin2dec(fliplr(dec2bin(j,bits))); end
        s=s(r(r<n)+1);
    otherwise
        error('pole:order','Unknown pole order %s.',order);
end
end
