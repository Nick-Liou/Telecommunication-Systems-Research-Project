
% Test the Gray Mapping 
% Calculate the Gp for different order constellations

% The name of the file we save the values 
nameOfTheFile = 'GpValues.mat' ;


distance = 2;
maxN = 10;
%plot = false; % If the plot is false then it loads (instead of cuclulating) some Gp values 
plot = true;

binaryGrayCodePrint = true ; 

if  exist( nameOfTheFile , 'file') == 2 && ~ plot 
        %This means the file exists :) 
        load(nameOfTheFile);
        fprintf('\nLoaded the file named : %s with pre-calculated values for the Gp \n' , nameOfTheFile);
        %nameOfTheFile % just to print the name of the file 
else 
        GpValues = zeros ( maxN, 2 );
end

ExternalSymbols = zeros ( maxN , 7);

for n = 2:maxN
        
        if size(GpValues,1)>= n  
                %Checks if it is inside the bounds of the GpValues
                if GpValues(n,1) ~= 0
                        %Checks if it was already calculated
                        fprintf('Skiped the calculation of the order n = %d, it was already calculated\n' , n);
                        continue
                end
        end
        
        tic
       
        m=2^n;
        
       
        
        if rem(n,2)==0
                a = 2^(n/2);
                maxRows = a+a/2-1;
                maxColumns = a;
                %GpValues(n,?)= ( 4/3*(a-2)^2+147/30*(a-2)+14/3 ) /m ;    % old delete 
                
                GpValues(n,2)= ( 40*(a-2)^2+147*(a-2)+140 ) /(m*30) ; % New better :)
        elseif n==3                
                maxRows = 4;
                maxColumns = 3;                
        else
                a = 2^((n-3)/2);
                maxRows = 3*a+3*a/2-1;
                maxColumns = 3*a;
        end
        
        
        [SymbolCoordinates,SymbolData,constellationVector , ~ ,constPower] = RegularHQAM(n,distance) ;
        
        
        if plot 
                %Ploting Gray Mapping
                scatterplot(constellationVector,[],[],'r*'); 
                grid
                drawnow
                hold on
                offset = 0;
                if ( rem(n,2) == 0)
                        for j=0:a-1
                                if(j~=0 && rem(j,2)==0)
                                        offset = offset+1;    
                                end
                                for k=0:a-1
                                        if binaryGrayCodePrint
                                                % Binary                                        
                                                text(real(SymbolCoordinates(k-offset+1+a/2-1,j+1)), imag(SymbolCoordinates(k-offset+1+a/2-1,j+1))+0.2, string(dec2bin(SymbolData(k-offset+1+a/2-1,j+1),n))); 
                                        else
                                                % Decimal
                                                text(real(SymbolCoordinates(k-offset+1+a/2-1,j+1)), imag(SymbolCoordinates(k-offset+1+a/2-1,j+1))+0.2, string(SymbolData(k-offset+1+a/2-1,j+1))); 
                                        end
                                end
                        end
                        
                elseif n==3
                        offset = 0;  
                        for j=0:2
                                if(j~=0 && rem(j,2)==0)
                                        offset = offset+1;    
                                end
                                for k=0:2
                                        if ( j==1 && k == 2) 
                                                continue
                                        end        
                                        if binaryGrayCodePrint
                                                % Binary
                                                text(real(SymbolCoordinates( k-offset+1+2-1,j+1)), imag(SymbolCoordinates(k-offset+1+2-1,j+1))+0.2, string(dec2bin(SymbolData( k-offset+1+2-1,j+1),n))); 
                                        else
                                                % Decimal
                                                text(real(SymbolCoordinates( k-offset+1+2-1,j+1)), imag(SymbolCoordinates(k-offset+1+2-1,j+1))+0.2, string(SymbolData( k-offset+1+2-1,j+1))); 
                                        end
                                               
                                end
                        end
                                                
                else
                        for j=0:3*a-1
                                if(j~=0 && rem(j,2)==0)
                                        offset = offset+1;    
                                end    

                                for k=0:3*a-1
                                        if ( ( j < a/2 || j > 3*a-1-a/2 ) && ( k < a/2 || k > 3*a-1-a/2 )  )
                                               continue 
                                        end
                                        
                                        if binaryGrayCodePrint                                                
                                                % Binary
                                                text(real(SymbolCoordinates( k-offset+1+3*a/2-1   , j+1)), imag(SymbolCoordinates( k-offset+1+3*a/2-1   , j+1))+0.2, string(dec2bin(SymbolData( k-offset+1+3*a/2-1   , j+1),n)));
                                        else
                                                % Decimal
                                                text(real(SymbolCoordinates( k-offset+1+3*a/2-1   , j+1)), imag(SymbolCoordinates( k-offset+1+3*a/2-1   , j+1))+0.2, string(SymbolData( k-offset+1+3*a/2-1   , j+1))); 
                                        end
                                end
                        end
                end
                
        end
        
        
        Gp=0;
        GpAccurate=0;

        for row=1:maxRows
                for columns=1:maxColumns
                        Gps = 0 ;
                        Ns = 0 ; 
                        if(SymbolData(row,columns) == -1 )
                                continue
                        end
                        S = SymbolData(row,columns); 
                        if row >1
                                %We go left
                                if SymbolData(row-1,columns) ~= -1
                                        Gps = Gps + SetBits(  bitxor(S ,SymbolData(row-1,columns)   ) );
                                        Ns = Ns + 1 ; 
                                end

                        end
                        if row < maxRows
                                %We go right
                                if SymbolData(row+1,columns) ~= -1
                                        Gps = Gps + SetBits(  bitxor(S ,SymbolData(row+1,columns)   ) );
                                        Ns = Ns + 1 ; 
                                end
                        end
                        if columns > 1
                                %We go left and down (diagonaly) 
                                if SymbolData(row,columns-1) ~= -1
                                        Gps = Gps + SetBits(  bitxor(S ,SymbolData(row,columns-1)   ) );
                                        Ns = Ns + 1 ; 
                                end
                        end
                        if columns < maxColumns
                                %We go up and right (diagonaly) 
                                if SymbolData(row,columns+1) ~= -1
                                        Gps = Gps + SetBits(  bitxor(S ,SymbolData(row,columns+1)   ) );
                                        Ns = Ns + 1 ; 
                                end
                        end

                        if row < maxRows && columns > 1
                                %We go down and right (diagonaly) 
                                if SymbolData(row+1,columns-1) ~= -1
                                        Gps = Gps + SetBits(  bitxor(S ,SymbolData(row+1,columns-1)   ) );
                                        Ns = Ns + 1 ; 
                                end
                        end
                        if row >1 && columns < maxColumns
                                %We go up and left (diagonaly) 
                                if SymbolData(row-1,columns+1) ~= -1
                                        Gps = Gps + SetBits(  bitxor(S ,SymbolData(row-1,columns+1)   ) );
                                        Ns = Ns + 1 ; 
                                end
                        end
                        
                        
                        GpsAccurate = Gps * (60 / Ns) ; 
                        GpAccurate = GpAccurate + GpsAccurate ;
                        
                        if Ns == 1                                 
                                ExternalSymbols(n,1) = ExternalSymbols(n,1) + 1 ;                       
                        elseif Ns == 2                                
                                ExternalSymbols(n,2) = ExternalSymbols(n,2) + 1 ;                       
                        elseif Ns == 3                                
                                ExternalSymbols(n,3) = ExternalSymbols(n,3) + 1 ;                       
                        elseif Ns == 4                                
                                ExternalSymbols(n,4) = ExternalSymbols(n,4) + 1 ;                       
                        elseif Ns == 5                                
                                ExternalSymbols(n,5) = ExternalSymbols(n,5) + 1 ;                       
                        elseif Ns == 6                                
                                ExternalSymbols(n,6) = ExternalSymbols(n,6) + 1 ;
                        end
                        % Not that accurate due to division errors
                        %Gps = Gps/Ns;
                        %Gp = Gp + Gps; 


                end
        end
        
        % Not that accurate due to division errors
        %Gp = Gp/m;
        %GpValues(n,?)= Gp;
        
        GpAccurate = GpAccurate / 60 / m ;
        GpValues(n,1)= GpAccurate;
        
        
        toc
        
        
end

ExternalSymbols(:,7) = ExternalSymbols(:,1) + ExternalSymbols(:,2) + ExternalSymbols(:,3) + ExternalSymbols(:,4) + ExternalSymbols(:,5) ;
GpValues %#ok<NOPTS>
ExternalSymbols;        


save(nameOfTheFile,'GpValues');


