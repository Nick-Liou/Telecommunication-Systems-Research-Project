function [SER , BER , GrayCode] = DetectionLinearSlowTesting(constellationVector , constellationGrayCodeVector , symbolsWithNoise , grayCodeOfSymbolSent  , m , d  )
        
        % This variation has the GrayCode enabled to run the DetectionTesting algorithm 
        
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
        
        % This -2 is used to work with the DetectionTesting 
        if ( grayCodeOfSymbolSent == -2 )
          	GrayCode = constellationGrayCodeVector(index) ;                         
                SER = 0 ;
            	BER = 0 ;
        else

                if grayCodeOfSymbolSent == constellationGrayCodeVector(index)
                        SER = 0 ;
                        BER = 0 ;

                        GrayCode = constellationGrayCodeVector(index) ;                         
                else
                        SER = 1 ; 
                        BER = SetBits( bitxor( grayCodeOfSymbolSent , constellationGrayCodeVector(index) )) ;
                        GrayCode = constellationGrayCodeVector(index) ;                         
                end
        end

end

