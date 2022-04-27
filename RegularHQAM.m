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
                
                
                GrayCodex = LinearGrayCode((n+1)/2); %3a x
                GrayCodey = LinearGrayCode((n-1)/2); %2a y
                
                %{
                GrayCode2D = zeros( 2^((n+1)/2) , 2^((n-1)/2) );
                for i=1:2^((n+1)/2)
                        for j= 1: 2^((n-1)/2)
                                GrayCode2D(i,j) =  bitshift( GrayCodex(i) , (n-1)/2  ) + GrayCodey(j);
                        end
                end
                GrayCode2D
                %}
                
                i = 1;                                                          	%matrix index
                offset = 0;  
                h_o = 3*a*distance/2 - distance/4;                                      %horizontal offset    +++++++
                v_o = (3*a/2 - 1)*sqrt(3)*distance/2 + sqrt(3)*distance/4;      	%vertical offset          +++++++

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

                                SymbolCoordinates(      k-offset+1+3*a/2-1        ,       j+1) = e(1) + 1i*e(2);
                                SymbolCoordinates2(     k-offset+1+3*a/2-1        ,       j+1) = inphase(i) +1i*quadr(i); 
                                
                                
                                %Gray code 
                                SymbolData(             k-offset+1+3*a/2-1        ,       j+1) = i+100000;       %starts from 1 ***Πρεπει να γινει Gray*** It should be uselles now !!! 
                                
                                if ( j >= a/2 && j <= 3*a-1-a/2 )
                                        
                                        SymbolData( k-offset+1+3*a/2-1 , j+1) =   GrayCodey(j+1-a/2) + bitshift( GrayCodex(k+1 + a/2) , (n-1)/2  )  ;  
                                        %Delete this error throw when we
                                        %make sure it works correct
                                        if bitand(  GrayCodey(j+1-a/2) , bitshift( GrayCodex(k+1 + a/2) , (n-1)/2  ) ) ~= 0
                                                errID = 'LinearGrayCode:TooBigInput';
                                                msg = 'Gray mapping error wrong bitshift !!!' ;
                                                exceptionWhenNisTooBig = MException(errID,msg);
                                                throw(exceptionWhenNisTooBig)
                                        end
                                else
                                        %The squares that are 'moved' from
                                        %the left and right of the 3a*2a
                                        % To swap we need the funtion f(i)= -i + x1 + x2 , where i=[x1,x2]  and f(i)=[x2,x1]
                                        %GrayCode2D
                                	if     ( j<a/2 &&  k < 3*a/2 )
                                                %Down left rectangle
                                                %SymbolData( k-offset+1+3*a/2-1   , j+1);
                                                
                                                %We separate it to 2 squares
                                                if ( k < 3*a/2 -a/2 )
                                                        %Left square Flip over the x axis
                                                        
                                                        %SymbolData( k-offset+1+3*a/2-1   , j+1) = GrayCode2D ( k-a/2+1 , j+1  ) ; % using the Gray2D before fliping over x
                                                        SymbolData( k-offset+1+3*a/2-1   , j+1) = bitshift( GrayCodex(  k-a/2+1 ) ,  (n-1)/2  ) + GrayCodey(a/2+1 - (j+1) ); 
                                                        
                                                else
                                                        %Right square over the y axis
                                                        %SymbolData( k-offset+1+3*a/2-1   , j+1) = GrayCode2D ( k-a+1 , j+1+a/2  ) ; % using the Gray2D before fliping over y
                                                        SymbolData( k-offset+1+3*a/2-1   , j+1) = bitshift( GrayCodex( a/2+1 -( k-a+1) ) ,  (n-1)/2  ) + GrayCodey( j+1+a/2  );
                                                end

                                        elseif ( j<a/2 && k >= 3*a/2 )
                                                %Down right rectangle
                                                %SymbolData( k-offset+1+3*a/2-1   , j+1);
                                                 if ( k<  2*a  )
                                                         %Left square over the y axis
                                                         %SymbolData( k-offset+1+3*a/2-1   , j+1) = GrayCode2D ( (k+2*a+1) , (j+a/2+1)  ) ; % using the Gray2D before fliping over y
                                                         %SymbolData( k-offset+1+3*a/2-1   , j+1) = GrayCode2D ( 15/2*a+1 - (k+2*a+1) , (j+a/2+1)  ) ; % using the Gray2D 
                                                         SymbolData( k-offset+1+3*a/2-1   , j+1) = bitshift( GrayCodex( 15/2*a+1 - (k+2*a+1) ) ,  (n-1)/2  ) + GrayCodey( (j+a/2+1)  );
                                                 else
                                                         %Right square  over the x axis
                                                         %SymbolData( k-offset+1+3*a/2-1   , j+1)  = GrayCode2D ( (k+a/2+1)+a , (j+a/2+1)-a/2   );  % using the Gray2D before fliping over x
                                                         %SymbolData( k-offset+1+3*a/2-1   , j+1)  = GrayCode2D ( (k+a/2+1)+a , 1+a/2 - ((j+a/2+1)-a/2)   );  % using the Gray2D 
                                                         SymbolData( k-offset+1+3*a/2-1   , j+1) = bitshift( GrayCodex( (k+a/2+1)+a ) ,  (n-1)/2  ) + GrayCodey(  1+a/2 - ((j+a/2+1)-a/2)  );
                                                         
                                                         
                                                 end
                                                
                                        elseif ( j >= 3*a/2 && k >= 3*a/2 )
                                                %Up right rectangle
                                                %SymbolData( k-offset+1+3*a/2-1   , j+1);
                                                
                                                %We separate it to 2 squares
                                                if ( k<  2*a  ) 
                                                        %Left square over the y axis
                                                        %SymbolData( k-offset+1+3*a/2-1   , j+1) = GrayCode2D ( (k+2*a+1) , (j+1-3/2*a)  ) ; % using the Gray2D before fliping over y
                                                        %SymbolData( k-offset+1+3*a/2-1   , j+1) = GrayCode2D ( 15/2*a+1 -(k+2*a+1) , (j+1-3/2*a)  ) ; % using the Gray2D
                                                        SymbolData( k-offset+1+3*a/2-1   , j+1) = bitshift( GrayCodex( 15/2*a+1 -(k+2*a+1) ) ,  (n-1)/2  ) + GrayCodey( (j+1-3/2*a) );
                                                else
                                                        %Right square  over the x axis
                                                        %SymbolData( k-offset+1+3*a/2-1   , j+1)  = GrayCode2D ( (k+a/2+1)+a , (j+1-3/2*a)+a/2  );  % using the Gray2D before fliping over x
                                                        %SymbolData( k-offset+1+3*a/2-1   , j+1)  = GrayCode2D ( (k+a/2+1)+a , 7/2*a+1-(j+1-3/2*a+a/2)  );  % using the Gray2D 
                                                        SymbolData( k-offset+1+3*a/2-1   , j+1) = bitshift( GrayCodex( (k+a/2+1)+a ) ,  (n-1)/2  ) + GrayCodey( 7/2*a+1-(j+1-3/2*a+a/2) );
                                                        
                                                end
                                                
                                        elseif ( j >= 3*a/2 && k < 3*a/2 )
                                                %Up left rectangle
                                                %SymbolData( k-offset+1+3*a/2-1   , j+1);
                                                
                                                %We separate it to 2 squares
                                                if ( k < 3*a/2 -a/2 )
                                                        %Left square over the x axis
                                                        %SymbolData( k-offset+1+3*a/2-1   , j+1) = GrayCode2D ( k-a/2+1 , (j+1-a)  ) ; % using the Gray2D before fliping over x
                                                        %SymbolData( k-offset+1+3*a/2-1   , j+1) = GrayCode2D ( k-a/2+1 , 7/2*a+1- (j+1-a)  ) ; % using the Gray2D 
                                                        SymbolData( k-offset+1+3*a/2-1   , j+1) = bitshift( GrayCodex(  k-a/2+1 ) ,  (n-1)/2  ) + GrayCodey( 7/2*a+1 -(j+1-a) );  
                                                else
                                                        %Right square  over the y axis
                                                        %SymbolData( k-offset+1+3*a/2-1   , j+1) = GrayCode2D ( k-a+1 , (j+1-3/2*a)  ) ; % using the Gray2D before fliping over y
                                                        SymbolData( k-offset+1+3*a/2-1   , j+1) = bitshift( GrayCodex( a/2+1 - (k-a+1) ) ,  (n-1)/2  ) + GrayCodey( (j+1-3/2*a) );
                                                end
                                	end


                                end
                                
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

