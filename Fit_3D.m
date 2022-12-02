
k_50 = fitresult1.k;
k_100 = fitresult3.k;
k_200 = fitresult2.k;

a_50 = fitresult1.a;
a_100 = fitresult3.a;
a_200 = fitresult2.a;

b_50 = fitresult1.b;
b_100 = fitresult3.b;
b_200 = fitresult2.b;

f_range = 50:1:300;
B_range = [0.05 0.1 0.2];

k = [k_50 k_100 k_200];
a = [a_50 a_100 a_200];
b = [b_50 b_100 b_200];

for B = 1:numel(B_range)
   for f = 1:numel(f_range)
        P(B,f) = k(B)*f_range(f)^a(B)*B_range(B)^b(B);  
   end
end
%P = k*x^a*0.1^b