function [SER , BER ] = DetectionLinearSlow(constellationVector , constellationGrayCodeVector , symbolsWithNoise , grayCodeOfSymbolSent  , m , d  )

        minDist = abs( constellationVector(1) - symbolsWithNoise);
        index=1;
        
        for i=2:m
                if( abs( constellationVector(i) - symbolsWithNoise) < minDist )
                        minDist = abs( constellationVector(i) - symbolsWithNoise);
                        index = i ;
                        
                        if (minDist < d/2)
                                break;
                        end
                        
                end 
        end
        
        if grayCodeOfSymbolSent == constellationGrayCodeVector(index)
                SER = 0 ;
                BER = 0 ;
                %GrayCode = constellationGrayCodeVector(index) ;
        else
                SER = 1 ; 
                BER = SetBits( bitxor( grayCodeOfSymbolSent , constellationGrayCodeVector(index) )) ;
                %GrayCode = constellationGrayCodeVector(index) ;
        end
       

end

