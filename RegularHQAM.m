function [SymbolCoordinates, SymbolData, constellationVector , constellationGrayCodeVector ,constPower] = RegularHQAM(n,distance)
        %[SymbolCoordinatesBeforeOffset,SymbolCoordinates2,SymbolCoordinates2Transpose,SymbolData,constellationVector , constellationGrayCodeVector ,constPower]
        
        %distance the minimum distance between 2 symbols
        %n order of the constellation
        
        m = 2^n;               	%number of symbols
        constellationGrayCodeVector = 1:m; 
        
        % A usefull constant
        if rem(n,2)==0
                a = 2^(n/2);                    
        else
                a = 2^((n-3)/2);
        end
        
        %SymbolCoordinatesBeforeOffset = [];    	%without offset                     
        SymbolCoordinates = [];    	%with offset - final coordinates     
        SymbolData = [];                %decimal number of symbol

        constPower = 0;         % The Mean Energy of the constellation (sum of squares)
        
        % Create basis vectors
        Q = [distance 0];
        R = [distance/2 sqrt(3)*distance/2];
       
        %preallocate the matrices for the symbols
        inphase = zeros(1,m);                           %x
        quadr = zeros(1,m);                             %y

        % Generate HQAM
        if rem(n,2)==0 && n>=2
                %for n even 
                %SymbolCoordinatesBeforeOffset = (zeros(a+a/2-1,a));          %without offset
                SymbolCoordinates = (zeros(a+a/2-1,a));     	%with offset - final coordinates
                SymbolData = (zeros(a+a/2-1,a))-1;            	%decimal number of symbol
              
                
                GrayCode = LinearGrayCode(n/2);
                if (isa( GrayCode(1) , 'int16') && n>=16) || (isa( GrayCode(1) , 'int32') && n>=32)
                        %Error with the bitshifting 
                        errID = 'RegularHQAM:TooSmallTypeForGrayStorage';
                        msg = 'In RegularHQAM at Gray mapping the type for storage is too small for that grey code, n='+ string(n)+'data type=' + class(GrayCode(1)) ;
                        exceptionWhenNisTooBig = MException(errID,msg);
                        throw(exceptionWhenNisTooBig)
                        
                end
                

                i = 1;        	%matrix index
                offset = 0;  
                h_o = a*distance/2 - distance/4;                               	%horizontal offset
                v_o = (a/2 - 1)*sqrt(3)*distance/2 + sqrt(3)*distance/4;     	%vertical offset
                for j=0:a-1
                        if(j~=0 && rem(j,2)==0)
                                offset = offset+1;    
                        end
                        for k=0:a-1
                                q = (k-offset)*Q;                               %x
                                r = j*R;                                        %y
                                e = q+r;                                        %symbol
                                inphase(i) = e(1) - h_o;
                                quadr(i) = e(2) - v_o;
                                SymbolData(k-offset+1+a/2-1,j+1) =  GrayCode(k+1) + bitshift( GrayCode(j+1) , n/2 );                           
                                %SymbolCoordinatesBeforeOffset(k-offset+1+a/2-1,j+1) = e(1) + 1i*e(2);
                                SymbolCoordinates(k-offset+1+a/2-1,j+1) = inphase(i) +1i*quadr(i);
                                constellationGrayCodeVector(i) = SymbolData(k-offset+1+a/2-1,j+1) ; % GrayCode 
                                i=i+1; 
                        end
                end
        elseif rem(n,2)~=0 && n>=5    
                %for n odd
                %SymbolCoordinatesBeforeOffset = (zeros(3*a+3*a/2-1 ,    3*a));  	%without offset
                SymbolCoordinates = (zeros(3*a+3*a/2-1 ,    3*a));  	%with offset - final coordinates
                SymbolData = (zeros(3*a+3*a/2-1 ,    3*a))-1;         	%decimal number of symbol
                
                
                GrayCodex = LinearGrayCode((n+1)/2); %3a x
                GrayCodey = LinearGrayCode((n-1)/2); %2a y
                
                if (isa( GrayCodex(1) , 'int16') && n>=16) || (isa( GrayCodex(1) , 'int32') && n>=32)
                        %Error with the bitshifting 
                        errID = 'RegularHQAM:TooSmallTypeForGrayStorage';
                        msg = 'In RegularHQAM at Gray mapping the type for storage is too small for that grey code, n='+ string(n)+'data type=' + class(GrayCodex(1)) ;
                        exceptionWhenNisTooBig = MException(errID,msg);
                        throw(exceptionWhenNisTooBig)
                        
                end
                
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
                h_o = 3*a*distance/2 - distance/4;                                      %horizontal offset              
                v_o = (3*a/2 - 1)*sqrt(3)*distance/2 + sqrt(3)*distance/4;      	%vertical offset                

                for j=0:3*a-1
                        if(j~=0 && rem(j,2)==0)
                                offset = offset + 1;    
                        end    

                        for k=0:3*a-1
                                if ( ( j < a/2 || j > 3*a-1-a/2 ) && ( k < a/2 || k > 3*a-1-a/2 )  )
                                	% Those are the egdes that don't have symbols
                                	continue
                                end
                                q = (k-offset)*Q;                               %x
                                r = j*R;                                        %y
                                e = q+r;                                        %symbol
                                inphase(i) = e(1) - h_o;
                                quadr(i) = e(2) - v_o;  

                                %SymbolCoordinatesBeforeOffset(      k-offset+1+3*a/2-1        ,       j+1) = e(1) + 1i*e(2);
                                SymbolCoordinates(     k-offset+1+3*a/2-1        ,       j+1) = inphase(i) +1i*quadr(i); 
                                
                                
                                %Gray code 
                                if ( j >= a/2 && j <= 3*a-1-a/2 )
                                        % Those are the easy ones, the center rectangle 3a*2a
                                        SymbolData( k-offset+1+3*a/2-1 , j+1) =   GrayCodey(j+1-a/2) + bitshift( GrayCodex(k+1 + a/2) , (n-1)/2  )  ;  
                                        
                                else
                                        %The squares that are 'moved' from the left and right of the 3a*2a
                                        % To swap we need the funtion f(i)= -i + x1 + x2 , where i=[x1,x2]  and f(i)=[x2,x1]
                                        %GrayCode2D
                                	if     ( j<a/2 &&  k < 3*a/2 )
                                                %Down left rectangle
                                                
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
                                                 if ( k<  2*a  )
                                                         %Left square over the y axis
                                                         %SymbolData( k-offset+1+3*a/2-1   , j+1) = GrayCode2D ( (k+2*a+1) , (j+a/2+1)  ) ; % using the Gray2D before fliping over y
                                                         %SymbolData( k-offset+1+3*a/2-1   , j+1) = GrayCode2D ( 15/2*a+1 - (k+2*a+1) , (j+a/2+1)  ) ; % using the Gray2D 
                                                         SymbolData( k-offset+1+3*a/2-1   , j+1) = bitshift( GrayCodex( 15/2*a+1 - (k+2*a+1) ) ,  (n-1)/2  ) + GrayCodey( (j+a/2+1)  );
                                                 else
                                                         %Right square  over the x axis
                                                         %SymbolData( k-offset+1+3*a/2-1   , j+1)  = GrayCode2D ( (k+a/2+1)+a , (j+a/2+1)-a/2   );  % using the Gray2D before fliping over x
                                                         %SymbolData( k-offset+1+3*a/2-1   , j+1)  = GrayCode2D ( (k+a/2+1)+a , 1+a/2 - ((j+a/2+1)-a/2)   );  % using the Gray2D 
                                                         %SymbolData( k-offset+1+3*a/2-1 , j+1) = bitshift( GrayCodex( (k+a/2+1)+a ) , (n-1)/2 ) + GrayCodey( 1+a/2 - ((j+a/2+1)-a/2) );  % The not simple one 
                                                         SymbolData( k-offset+1+3*a/2-1 , j+1) = bitshift( GrayCodex( k+3/2*a+1 ) , (n-1)/2 ) + GrayCodey( 1+a/2 - (j+1) );
                                                         
                                                 end
                                                
                                        elseif ( j >= 3*a/2 && k >= 3*a/2 )
                                                %Up right rectangle
                                                
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
                                                        SymbolData( k-offset+1+3*a/2-1   , j+1) = bitshift( GrayCodex( (k+a/2+1)+a ) ,  (n-1)/2  ) + GrayCodey( 7/2*a+1-(j+1-a) );
                                                        
                                                end
                                                
                                        elseif ( j >= 3*a/2 && k < 3*a/2 )
                                                %Up left rectangle
                                                
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
                                constellationGrayCodeVector(i) =  SymbolData( k-offset+1+3*a/2-1   , j+1) ;
                                i=i+1; 
                        end
                end
        elseif n==3
                
                
                %SymbolCoordinatesBeforeOffset = zeros(4,3);          %without offset
                SymbolCoordinates = zeros(4,3);       	%with offset - final coordinates
                SymbolData = (zeros(4,3))-1  ;      	%decimal number of symbol
              
                % Gray code hard coded 
                SymbolData(1,3) = 2;
                SymbolData(2,1) = 7;
                SymbolData(3,1) = 5;
                SymbolData(2,2) = 6;
                SymbolData(3,2) = 0;
                SymbolData(2,3) = 3;
                SymbolData(3,3) = 1;
                SymbolData(4,1) = 4;
                
                i = 1;                       	%matrix index
                offset = 0;  
                h_o = distance;            	%horizontal offset
                v_o = sqrt(3)*distance/2;  	%vertical offset
                for j=0:2
                        if(j~=0 && rem(j,2)==0)
                                offset = offset+1;    
                        end
                        for k=0:2
                                if ( j==1 && k == 2) 
                                        continue
                                end
                                q = (k-offset)*Q;                               %x
                                r = j*R;                                        %y
                                e = q+r;                                        %symbol
                                inphase(i) = e(1) - h_o;
                                quadr(i) = e(2) - v_o;                              
                                %SymbolCoordinatesBeforeOffset(k-offset+1+2-1,j+1) = e(1) + 1i*e(2);
                                SymbolCoordinates(k-offset+1+2-1,j+1) = inphase(i) +1i*quadr(i);
                                constellationGrayCodeVector(i) = SymbolData(k-offset+1+2-1,j+1);
                                i=i+1; 
                        end
                end

        end
        

        inphase = inphase(:);
        quadr = quadr(:);
        constellationVector = inphase + 1i*quadr;

        
        %Calculates the Constellation Power
        for i=1:m
                constPower = constPower + abs(constellationVector(i))^2; 
                %constPower =  real(constellationVector(i))^2 + imag(constellationVector(i))^2;
        end
        constPower=constPower/m;
        
        
end

