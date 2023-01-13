function [H]=CHxy_d(x,y)
nsamples=size(x,1);
Hxy=0;
%Sx£¬Sy
z=[x,y];
z2=sortrows(z,1);
[d,~,~]=unique(z2,'rows');
lend=size(d,1);
lenz=size(z,1);
temp=zeros(lend,1);
for i = 1:lend
    for j=1:lenz
        if d(i,:)==z2(j,:)
            temp(i)=temp(i)+1;
        end
    end
end
%H(X,Y)
for i=1:lend
Hxy=Hxy+(temp(i)/nsamples)*(log(temp(i)/nsamples));
end
H=-Hxy;
