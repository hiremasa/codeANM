%function ind = checkDSeparate(fIdx, condSetIdx, label, d, alpha,width,semi)%%%%%%%%%%%%%%%%%%%%%20130814
function ind = D_Separate(fIdx, condSetIdx, label, d, p_vaule)
%fIdx, condSetIdx下标
[nS, nF]=size(d);
nLS=length(label);
%%%%%%%%%%%%%%%%%%%%%20130814
%if intersect(fIdx,semi)
   % ind=false;%有关
   % return
%end
%%%%%%%%%%%%%%%%%%%%%20130814
%ind=false;
%%%%%%%%%% shrink the condition set using independence（利用了算法3）%%%%%%%%%%
%temp=false(length(condSetIdx),1);
%for i=1:length(condSetIdx)
  %  Ixy2 =  CMI(d(1:nS, fIdx), d(1:nS, condSetIdx(i)));
   % if Ixy2<p_vaule
   %     ind=true;
   % else ind=false;
   % end
   % temp(i)=ind;
%end;

%condSetIdx=condSetIdx(find(temp==false));%找到有关的
%%%%%%%%%% shrink the condition set using cond independenc（利用了算法3）%%%%%%%
temp=false(length(condSetIdx),1);
for i=1:length(condSetIdx)
    for j=1:length(condSetIdx)
        
        if temp(i)==false && i~=j
            Ixy3  = CMI(d(1:nS, fIdx), d(1:nS, condSetIdx(i)),d(1:nS, condSetIdx(j)));
            if Ixy3<p_vaule
                ind=true;
            else
                ind=false;
            end
            temp(i)=ind;
        end
    end
end;
condSetIdx=condSetIdx(find(temp==false));

if isempty(condSetIdx)
    Ixy5 =CMI(d(1:nLS, fIdx), label);
    if Ixy5<p_vaule
        ind=true;
    else
        ind=false;
    end
    return
else
    %%%%%%%%%% enumerate all possible subsets of condSetIdx（遍历pc节点的所有子集）%%%%%%%%
    len=length(condSetIdx);
    for i=1:2^len
        B=dec2bin(i,len);
        condSubSet=condSetIdx(find(B=='1'));
        if length(condSetIdx)<=4
            Ixy4=CMI(d(1:nLS, fIdx),label, d(1:nLS,condSubSet));
            if Ixy4<p_vaule
                ind=true;
                return ;
            end
        end
    end
end;
ind=false;
