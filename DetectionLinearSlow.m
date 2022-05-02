function [SER , BER ] = DetectionLinearSlow(constellationVector , constellationGrayCodeVector , symbolsWithNoise , grayCodeOfSymbolSent  , m , d  )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

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
        
        SER=0;
        BER = SetBits( bitxor( grayCodeOfSymbolSent , constellationGrayCodeVector(index) )) ;
        if (BER ~= 0 )
                SER=1;
        end
        
                


end

