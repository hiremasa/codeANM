 
function  final_result=perm_comb(x,y)
m=1:x;
n=[];
temp=combnk(m,y);
for k=1:size(temp,1)
    n=[n;perms(temp(k,:))];
end
aa=repmat(m(:),1,y);
final_result=[aa;n];
