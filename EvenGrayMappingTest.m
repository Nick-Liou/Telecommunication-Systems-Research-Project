
% Test the Gray Mapping 

% for n =16 we have a problem with the SymbolData there are NEGATIVE values

distance = 2;
maxN = 18;
GpValues = zeros ( maxN, 1 );
plot = false; 
%plot = true;

for n = 2:maxN
        m=2^n;
        
        if n ==3 
                continue % We have not made this yet
        end
        
        if rem(n,2)==0
                a = 2^(n/2);    
                % a+a/2-1,a
                maxRows = a+a/2-1;
                maxColumns = a;
        else
                a = 2^((n-3)/2);
                %3*a+3*a/2-1 ,    3*a
                maxRows = 3*a+3*a/2-1;
                maxColumns = 3*a;
        end
        
        [~,SymbolCoordinates2,~,SymbolData,refConst,~] = RegularHQAM(n,distance) ;
        
        if plot 
                %Ploting Gray Mapping
                scatterplot(refConst,[],[],'r*');
                grid
                drawnow
                hold on
                offset = 0;
                for j=0:a-1
                        if(j~=0 && rem(j,2)==0)
                                offset = offset+1;    
                        end
                        for k=0:a-1
                                %text(real(SymbolCoordinates2(k-offset+1+a/2-1,j+1)), imag(SymbolCoordinates2(k-offset+1+a/2-1,j+1))+0.2, string(dec2bin(SymbolData(k-offset+1+a/2-1,j+1),n))); 
                                text(real(SymbolCoordinates2(k-offset+1+a/2-1,j+1)), imag(SymbolCoordinates2(k-offset+1+a/2-1,j+1))+0.2, string(SymbolData(k-offset+1+a/2-1,j+1))); 
                        end
                end
        end
        
        
        Gp=0;

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
                        Gps = Gps/Ns;
                        Gp = Gp + Gps; 


                end
        end

        Gp = Gp/m;
        GpValues(n)= Gp;

        
        
end

GpValues

%{

%n=6;
m=2^n;
distance = 2;

if rem(n,2)==0
        a = 2^(n/2);    
        % a+a/2-1,a
        maxRows = a+a/2-1;
        maxColumns = a;
else
        a = 2^((n-3)/2);
        %3*a+3*a/2-1 ,    3*a
        maxRows = 3*a+3*a/2-1;
        maxColumns = 3*a;
end

%[SymbolCoordinates,SymbolCoordinates2,SymbolCoordinates2Transpose,SymbolData,refConst,constPower] = RegularHQAM(n,distance) ;
[~,SymbolCoordinates2,~,SymbolData,refConst,~] = RegularHQAM(n,distance) ;

SymbolData %#ok<NOPTS> % Printing The Gray Dec values



Gp=0;

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
                Gps = Gps/Ns;
                Gp = Gp + Gps; 


        end
end

Gp = Gp/m






%Ploting Gray Mapping
scatterplot(refConst,[],[],'r*');
grid
drawnow
hold on
offset = 0;
for j=0:a-1
        if(j~=0 && rem(j,2)==0)
                offset = offset+1;    
        end
        for k=0:a-1
                %text(real(SymbolCoordinates2(k-offset+1+a/2-1,j+1)), imag(SymbolCoordinates2(k-offset+1+a/2-1,j+1))+0.2, string(dec2bin(SymbolData(k-offset+1+a/2-1,j+1),n))); 
                text(real(SymbolCoordinates2(k-offset+1+a/2-1,j+1)), imag(SymbolCoordinates2(k-offset+1+a/2-1,j+1))+0.2, string(SymbolData(k-offset+1+a/2-1,j+1))); 
        end
end

%}