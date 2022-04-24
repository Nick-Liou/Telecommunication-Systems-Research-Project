function [SymbolCoordinates,SymbolCoordinates2,SymbolCoordinates2Transpose,SymbolData,refConst,constPower] = RegularHQAM(n,distance)
%UNTITLED6 Summary of this function goes here ++++++++
%   Detailed explanation goes here ++++++++++++++
        %inputs
        %distance the minimum distance between 2 symbols
        %n order of the constellation
        
        m = 2^n;               	%number of symbols
        
        
        % A usefull constant
        if rem(n,2)==0
                a = 2^(n/2);                    
        else
                a = 2^((n-3)/2);
        end


        SymbolCoordinates = [];    	%without offset
        SymbolCoordinates2 = [];    	%with offset - final coordinates
        SymbolData = [];                %decimal number of symbol

        constPower = 0;         % The Mean Energy of the constellation (sum of squares)
        
        % Do we need vectors ? or just complex numbers ?
        % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        % Create basis vectors
        Q = [distance 0];
        R = [distance/2 sqrt(3)*distance/2];
        %Q = Q.';
        %R = R.';

        %preallocate the matrices for the symbols
        inphase = zeros(1,m);                           %x
        quadr = zeros(1,m);                             %y

        % Generate HQAM
        if rem(n,2)==0 && n>=2
                SymbolCoordinates = (zeros(a+a/2-1,a));          %without offset
                SymbolCoordinates2 = (zeros(a+a/2-1,a));        %with offset - final coordinates
                SymbolData = (zeros(a+a/2-1,a))-1;                       %decimal number of symbol
              
                
                GrayCode = LinearGrayCode(n/2);

                i = 1;                                                                           %matrix index
                offset = 0;  
                h_o = a*distance/2 - distance/4;                                                  %horizontal offset
                v_o = (a/2 - 1)*sqrt(3)*distance/2 + sqrt(3)*distance/4;          %vertical offset
                for j=0:a-1
                        if(j~=0 && rem(j,2)==0)
                                offset = offset+1;    
                        end
                        for k=0:a-1
                                q = (k-offset)*Q;                       %x
                                r = j*R;                                        %y
                                e = q+r;                                        %symbol
                                inphase(i) = e(1) - h_o;
                                quadr(i) = e(2) - v_o;
                                SymbolData(k-offset+1+a/2-1,j+1) =          GrayCode(k+1) + bitshift( GrayCode(j+1) , n/2 )    ;                           
                                SymbolCoordinates(k-offset+1+a/2-1,j+1) = e(1) + 1i*e(2);
                                SymbolCoordinates2(k-offset+1+a/2-1,j+1) = inphase(i) +1i*quadr(i);
                                i=i+1; 
                        end
                end
        elseif rem(n,2)~=0 && n>=5    
                %for n odd
                SymbolCoordinates = (zeros(3*a+3*a/2-1 ,    3*a));     %without offset
                SymbolCoordinates2 = (zeros(3*a+3*a/2-1 ,    3*a));    %with offset - final coordinates
                SymbolData = (zeros(3*a+3*a/2-1 ,    3*a))-1;            %decimal number of symbol

                i = 1;                                                                           %matrix index
                offset = 0;  
                h_o = 3*a*distance/2 - distance/4;                                                  %horizontal offset    +++++++
                v_o = (3*a/2 - 1)*sqrt(3)*distance/2 + sqrt(3)*distance/4;          %vertical offset          +++++++

                for j=0:3*a-1
                        if(j~=0 && rem(j,2)==0)
                                offset = offset+1;    
                        end    

                        for k=0:3*a-1
                                if ( ( j < a/2 || j > 3*a-1-a/2 ) && ( k < a/2 || k > 3*a-1-a/2 )  )
                                       continue
                                end
                                q = (k-offset)*Q;                               %x
                                r = j*R;                                        %y
                                e = q+r;                                        %symbol
                                inphase(i) = e(1) - h_o;
                                quadr(i) = e(2) - v_o;

                                SymbolData(             k-offset+1+3*a/2-1        ,       j+1) = i;       %starts from 1 ***Πρεπει να γινει Gray***
                                SymbolCoordinates(      k-offset+1+3*a/2-1        ,       j+1) = e(1) + 1i*e(2);
                                SymbolCoordinates2(     k-offset+1+3*a/2-1        ,       j+1) = inphase(i) +1i*quadr(i); 
                                i=i+1; 
                        end
                end
        elseif n==3
                inphase = [0 distance distance/2 distance];
                quadr = [sqrt(3)*distance/4 sqrt(3)*distance/4 0 -sqrt(3)*distance/4];
                inphase = [inphase; -inphase];
                quadr = [quadr; -quadr];
                % NEEDs more work +++++++


        end



        SymbolCoordinates2Transpose = SymbolCoordinates2.';

        inphase = inphase(:);
        quadr = quadr(:);
        refConst = inphase + 1i*quadr;

        
        %Calculates the Constellation Power
        for i=1:m
                constPower = constPower + abs(refConst(i))^2; % real(refConst(i))^2 + imag(refConst(i))^2;
        end
        constPower=constPower/m;
        
end

