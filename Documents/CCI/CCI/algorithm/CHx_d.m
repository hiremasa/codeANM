function [H]=CHx_d(x)
nsamples=size(x,1);
Hx=0;
B=unique(x);
lenB=length(B);
for i=1:lenB
Hx=Hx+(length(find(x==B(i)))/nsamples)*(log((length(find(x==B(i)))/nsamples)));
end
Hx=-Hx;
H=Hx;
