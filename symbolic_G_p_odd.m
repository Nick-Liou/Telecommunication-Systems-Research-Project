

a = sym("a");

% A way to split the symbols
symbols = [
  (6*a-2)*(4*a-2) 
 8
 (4*a-2)
 (4*a-2)
 2*(2*a-3)
 (a-2)^2 * 8
 (a-2) * 2
 (a-2) * 2
 (a-2) * 8
 (a-2) * 3 * 2
 (a-2) * 3 * 2
 2     * 3 * 2 
 ( (4*a-3) + (4*a+1-3) ) * 2 
 2  * 2 
 2
 2
 2
 2
] ;

% Test that we counted all symbols
% it is correct since they should be 32*a^2
total_symbols = sum(symbols);
total_symbols = simplify(total_symbols)

% a weighted sum of the symbols to find the Gp_odd
G_p_s = [
4/3 *  (6*a-2)*(4*a-2) 
14/3 * 2
3/3 * (4*a-2)
7/5 * (4*a-2)
5/4 * 2*(2*a-3)
4/3 * (a-2)^2 * 8
3/3 * (a-2) * 2
7/5 * (a-2) * 2
5/4 * (a-2) * 8
4/3 * (a-2) * 3 * 2
4/3 * (a-2) * 3 * 2
5/4 * 2     * 3 * 2 
5/3 * ( (4*a-3) + (4*a+1-3) ) * 2 
4/3 * 2  * 2 
7/5 * 2
9/5 * 2
6/5 * 2
5/4 * 2
] ;


% Simplify the closed form of the formula
Gray_code_penatly = sum(G_p_s);
Gray_code_penatly = simplify(Gray_code_penatly)


Gp_odd = simplify(Gray_code_penatly / total_symbols)

% Works for n>=7  (n=3 and n=5 are special cases)
n = 17 ; 
a_num = 2 ^ ((n-5)/2) ; 
Gp_odd = subs(Gp_odd,a,a_num) ; 

fprintf("For n = %d , a = %d, " , n , a_num );
fprintf("Gp_odd = ");
disp(Gp_odd)
fprintf("Gp_odd = ");
disp(double(Gp_odd))




