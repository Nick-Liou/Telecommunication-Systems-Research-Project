function [codeDecimal] = LinearGrayCode(n)
% Returns a vector of 2^n elements with intigers gray coded
        if n>=32 
                errID = 'LinearGrayCode:TooBigInput';
                msg = 'LinearGrayCode input must be at most 31 and was given: '+ string(n) ;
                exceptionWhenNisTooBig = MException(errID,msg);
                throw(exceptionWhenNisTooBig)
        end
        
        codeDecimal = zeros(1,2^n, 'int32');
        
        for i = 0:2^n-1
                % i is the i'th number of the gray code table
                for j = 0:n-1
                        % j is the j'th bit of the gray code
                        a = idivide( mod( int32(i) , 2^(n-j) ) , 2^(n-j-1 ) );
                        b = mod( idivide( int32(i) , 2^(n-j) ) , int32(2) );
                        codeDecimal(i+1) = codeDecimal(i+1) + bitshift( abs(a-b) , n-j-1 );
                end
                
        end
        
end

