%Gey testing :)

%Max of n=15 because The Grey code is int16
errors =0;
for n=1:15 
        tic
        %n=12;
        GreyCode = LinearGreyCode(n);
        for i= 1:2^n-1
                if SetBits( bitxor( GreyCode(i),GreyCode(i+1) ))~= 1 
                        errors = errors +1;
                end
        end
        toc
end