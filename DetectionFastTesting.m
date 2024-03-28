function [SER , BER , GrayCode] = DetectionFastTesting(SymbolCoordinates , SymbolData , symbolWithNoise , grayCodeOfSymbolSent  , p , distance , a)

    % This variation has the GrayCode enabled to run the DetectionTesting algorithm 
       
    if p == 1
            % p=1                 
            h_o = a*distance/2 - distance/4;                               	%horizontal offset
            v_o = (a/2 - 1)*sqrt(3)*distance/2 + sqrt(3)*distance/4;     
    elseif p == 2
            %p=2
            h_o = 3*a*distance/2 - distance/4;                            	%horizontal offset              
            v_o = (3*a/2 - 1)*sqrt(3)*distance/2 + sqrt(3)*distance/4; 
    else
            %p==3
            h_o = distance;                                             	%horizontal offset
            v_o = sqrt(3)*distance/2;
    end
    
    
    
    s = distance / sqrt(3);
    x = real(symbolWithNoise) + h_o;
    y = imag(symbolWithNoise) + v_o;
    
    
    
    %y = 3/2 * s * b
    b = 2/3 * y / s;
    %x = sqrt(3) * s * ( b/2 + r)
    %x = - sqrt(3) * s * ( b/2 + g )
    r = (sqrt(3)/3 * x - y/3 ) / s;
    g = -(sqrt(3)/3 * x + y/3 ) / s;

    %r + b + g = 0
       
    %Rounding
    b_r = round(b);
    r_r = round(r);
    g_r = round(g);

    b_diff = abs(b_r - b);
    r_diff = abs(r_r - r);
    g_diff = abs(g_r - g);

    if (g_diff > r_diff) && (g_diff > b_diff)
        %g_r = -r_r - b_r;      %This is not used anywhere 
    elseif r_diff > b_diff
        r_r = -g_r - b_r;
    else
        b_r = -g_r-r_r;
    end  
    
    
    sz = size(SymbolCoordinates);
    if  p==1
            r_r_f = r_r + 1 + sz(1)-a  ;
            b_r_f = b_r + 1 ; 
            
    elseif p==2            
            r_r_f = r_r + 1 + 3/2*a -1 ;
            b_r_f = b_r + 1 ; 
           
    else % p==3        
           r_r_f = 2 + r_r ; 
           b_r_f = b_r + 1 ; 
          
    end
    
    
    
    % We have culculated the indexes of the symbol sent in the matrix 
    % Now we need to check if it is in the 2D matrix and if it is if there is a symbol in that cell 
    % If not we must find the closest symbol 
    
    try 
            % Tries to see if it is in the 2D matrix
            %symbolDetected = SymbolCoordinates( r_r_f , b_r_f);
            symbolDetectedGrayCode = SymbolData( r_r_f , b_r_f);

            if ( symbolDetectedGrayCode == grayCodeOfSymbolSent) 
                    %Symbol is in the correct hex :)
                    SER = 0 ;
                    BER = 0 ;
                    
                    GrayCode = symbolDetectedGrayCode ;
            else
                    % It is in the wrong hex
                    % Either that hex is 'valid' / of another symbol 
                    % OR it is from the 'blank' symbols so it has a gray code -1 ( and sybmol coordinates 0+1i*0 ) 
                    
                    % Try to evaluate the SetBits
                    % If the hex detected is in a blank hex ( so symbolDetectedGrayCode = -1 ) the bitxor will throw an error:
                    % Error using bitxor
                    % Double inputs must have integer values in the range of ASSUMEDTYPE.
                    % Because symbolDetectedGrayCode is -1

                    % If it isn't a 'blank' hex then it will compute correctly 
                    BER = SetBits( bitxor( symbolDetectedGrayCode , grayCodeOfSymbolSent)) ;
                    SER = 1 ; 

                    GrayCode = symbolDetectedGrayCode ;
                           
            end
    catch 
                % Here he have to check the neighbors                              
                % If the're isn't a non blank neighbor then call slow detect
                maxRows = sz(1) ;
                maxColumns = sz(2) ;
                row = r_r_f ; 
                columns = b_r_f;                 
                minDist = inf; 
                index_i = row ;
                index_j = columns ;
                
                if row >1        &&  row <=maxRows+1 && 1<=columns && columns <= maxColumns
                        %We go left
                        
                        if SymbolData(row-1,columns) ~= -1                                
                                 if( abs( SymbolCoordinates(row-1,columns) - symbolWithNoise) < minDist )                                
                                        minDist = abs( SymbolCoordinates(row-1,columns) - symbolWithNoise);
                                        index_i = row-1 ;
                                        index_j = columns ; % This is 'kinda' extra since it already has this value                                         
                                 end                                
                        end   
                end
                
                if row < maxRows         &&  row >=0 && 1<=columns && columns <= maxColumns
                        %We go right
                        if SymbolData(row+1,columns) ~= -1                                
                                if( abs( SymbolCoordinates(row+1,columns) - symbolWithNoise) < minDist )                                
                                        minDist = abs( SymbolCoordinates(row+1,columns) - symbolWithNoise);
                                        index_i = row+1 ;
                                        index_j = columns ; % This is 'kinda' extra since it already has this value                                         
                                end
                        end
                end
                
                if columns > 1          && row >=1  &&  row <=maxRows && columns<= maxColumns+1
                        %We go left and down (diagonaly) 
                        if SymbolData(row,columns-1) ~= -1                                
                                if( abs( SymbolCoordinates(row,columns-1) - symbolWithNoise) < minDist )                                
                                        minDist = abs( SymbolCoordinates(row,columns-1) - symbolWithNoise);
                                        index_i = row ;
                                        index_j = columns-1 ;                                        
                                end
                        end
                end
                
                if columns < maxColumns         && row >=1  &&  row <=maxRows && columns >= 0 
                        %We go up and right (diagonaly) 
                        if SymbolData(row,columns+1) ~= -1                                
                                if( abs( SymbolCoordinates(row,columns+1) - symbolWithNoise) < minDist )                                
                                        minDist = abs( SymbolCoordinates(row,columns+1) - symbolWithNoise);
                                        index_i = row ;
                                        index_j = columns+1 ;                                        
                                end
                        end
                end

                if row < maxRows && columns > 1         && row >=0 && columns <= maxColumns+1
                        %We go down and right (diagonaly) 
                        if SymbolData(row+1,columns-1) ~= -1                                
                                if( abs( SymbolCoordinates(row+1,columns-1) - symbolWithNoise) < minDist )                                
                                        minDist = abs( SymbolCoordinates(row+1,columns-1) - symbolWithNoise);
                                        index_i = row+1 ;
                                        index_j = columns-1 ;                                        
                                end
                        end
                end
                
                if row >1 && columns < maxColumns       && row <= maxRows+1 && columns >= 0
                        %We go up and left (diagonaly) 
                        if SymbolData(row-1,columns+1) ~= -1                                
                                if( abs( SymbolCoordinates(row-1,columns+1) - symbolWithNoise) < minDist )                                
                                        minDist = abs( SymbolCoordinates(row-1,columns+1) - symbolWithNoise);
                                        index_i = row-1 ;
                                        index_j = columns+1 ;                                        
                                end
                        end
                end
                        
                
                
                
                if minDist ~= inf 
                        if  grayCodeOfSymbolSent  == SymbolData( index_i , index_j )
                                SER = 0 ;
                                BER = 0 ;
                                GrayCode = grayCodeOfSymbolSent ;
                        else
                                SER = 1;
                                BER = SetBits(  bitxor(   grayCodeOfSymbolSent  ,SymbolData(   index_i , index_j       )   ) ) ;
                                GrayCode = SymbolData(index_i , index_j) ;
                        end
                else 
                        [SER , BER , GrayCode ] = DetectionLinearSlow2D(SymbolCoordinates , SymbolData , symbolWithNoise , grayCodeOfSymbolSent );
                        %[SER , BER , ~ ] = DetectionLinearSlow2D(SymbolCoordinates , SymbolData , symbolWithNoise , grayCodeOfSymbolSent );
                        
                end
                       
    end
    
    
    
end    
   

