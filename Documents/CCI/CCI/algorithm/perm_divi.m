 
function  final_result=perm_divi(n)
 a = rand;
 b = rand;
 c = rand;
 d = a+b+c;
 a = round(n*a/d);
 b = round(n*b/d);
 if a == 0
     a = 1;
 end
 if b == 0
     b = 1;
 end
  if a > 7
     a = 7;
  end
  if b > 7
     b = 7;
 end
 if a + b == n
     if a == max(a,b);
         a = a - 1;
     else
         b = b - 1;
     end
 end
 c = n - a - b;
 final_result = [a ,b ,c];
end
