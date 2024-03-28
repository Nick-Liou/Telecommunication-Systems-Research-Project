%Gray code testing :)


errors =0;
for n=1:20 
        tic
        
        GrayCode = LinearGrayCode(n);        
        for i= 1:2^n-1
                if SetBits( bitxor( GrayCode(i),GrayCode(i+1) ))~= 1 
                        errors = errors +1;
                end
        end
        
        toc
end

if errors ~= 0 
        fprintf('There was an error and the linear Gray code was created wrong \n' );
else        
        fprintf('The linear Gray code is correct \n' );
end 