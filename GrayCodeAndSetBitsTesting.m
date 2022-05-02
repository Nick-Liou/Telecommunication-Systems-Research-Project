%Gey testing :)

errors =0;
for n=1:15 
        tic
        %n=12;
        GrayCode = LinearGrayCode(n);
        for i= 1:2^n-1
                if SetBits( bitxor( GrayCode(i),GrayCode(i+1) ))~= 1 
                        errors = errors +1;
                end
        end
        toc
end