function [SER , BER , GrayCode ] = DetectionLinearSlow2D(SymbolCoordinates , SymbolData , symbolWithNoise , grayCodeOfSymbolSent )
        
        minDist = inf;
        index_i = -1;
        index_j = -1; 

        sizes = size(SymbolCoordinates);
        for i = 1:sizes(1)
                for j = 1:sizes(2)
                        if SymbolData(i,j) ~= -1 
                                if( abs( SymbolCoordinates(i,j) - symbolWithNoise) < minDist )
                                
                                        minDist = abs( SymbolCoordinates(i,j) - symbolWithNoise);
                                        index_i = i ;
                                        index_j = j ;
                                        %{
                                        % It is called in the DetectionFast when the symbol with noise is not close to any of the symbols in the constellation 
                                        % So it is not expected to be closer than d/2
                                        if (minDist < distance/2)
                                                break;
                                        end
                                        %}
                                
                                end
                        end
                end
        end

        BER = SetBits( bitxor( grayCodeOfSymbolSent , SymbolData(index_i,index_j) )) ;
        if BER ~= 0 
                SER = 1 ; 
        else
                SER = 0 ; 
        end
        GrayCode = SymbolData(index_i,index_j) ;
        
        
end

