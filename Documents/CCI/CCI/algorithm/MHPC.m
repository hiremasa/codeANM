
function  [trueCC] = MHPC(CC,data)
trueCC = [];
if isempty(data)
    fprintf('Warning, there is no remaining cause factor\n\n');
    return;
end
[ind1,~]  = Independent_test(CC, data(:,1), 0.05);
[ind2,~]  = Conditional_independent_test(CC, data(:,1), data(:,2), 0.05);
if ind1 == 1 && ind2 ==0
    trueCC = 1;
end
end
