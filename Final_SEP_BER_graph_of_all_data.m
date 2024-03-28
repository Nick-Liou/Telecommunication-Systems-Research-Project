
% Makes the final graph for SEP and BER

nameOfThePointsFile1 = "MultiRunsValuesForSEP&BER_Nick.mat" ;
nameOfThePointsFile2 = "MultiRunsValuesForSEP&BER_Katerina.mat" ;
nameOfThePointsFile3 = "MultiRunsValuesForSEP&BER_Nick_Old.mat" ;
nameOfThePointsFile4 = "MultiRunsValuesForSEP&BER_Extra.mat" ;


namesOfPointsFiles = [ nameOfThePointsFile1  nameOfThePointsFile2 nameOfThePointsFile3 nameOfThePointsFile4  ] ; 
dataLoaded = false; 

maxSNR = 0 ; 
maxN=0 ;

for file = namesOfPointsFiles
        if  exist( file , 'file') == 2 
                %This means the file exists :)
                load(file);
                fprintf('\nLoaded the file named : %s with last SEP and BER values ' , file);
                
                % Now we need to add the new data points to the old one 
                dim = size(MultiRunsValues) ; 
                %  n , snr , Number of Run ,  vector of data
                % vector of data :  = [ SEP BER NOS Time DetecAlgoID Who When ] 
                
                if ~dataLoaded 
                        % If it is the first time we load data
                	statsMatrix = zeros ( dim(1), dim(2),  4 );
                        statsMatrixLinearSlow =  zeros ( dim(1), dim(2),  4 );
                        statsMatrixFastDetection = zeros ( dim(1), dim(2),  4 );
                        
                        dataLoaded = true; 
                end
                

                %maxN = dim(1);
                if dim(2) > maxSNR
                        maxSNR =  dim(2); 
                end
                if dim(1) > maxN
                        maxN= dim(1);
                end
                
                for i= 1:dim(3)
                        statsMatrix =  statsMatrix + squeeze( MultiRunsValues ( : , : , i , 1:4) );
                        if max(max( MultiRunsValues ( : , : , i , 5 ) )) ==  10                              
                                statsMatrixLinearSlow =  statsMatrixLinearSlow + squeeze( MultiRunsValues ( : , : , i , 1:4) );
                        end
                        if max(max( MultiRunsValues ( : , : , i , 5 ) )) ==  20                             
                                statsMatrixFastDetection =  statsMatrixFastDetection + squeeze( MultiRunsValues ( : , : , i , 1:4) );
                        end
                end
        
        else
                fprintf('\nThe file named : %s does not exist in the current directory' , file);
        end
end


if dataLoaded 
        berFix = zeros(40,10);
        for i=1:40
                berFix(i,:) = 1:maxN;
        end
        
        snrvalues = 1:1:maxSNR;
        SEPvector = transpose ( statsMatrix(: , : ,1) ./ statsMatrix(: , : ,3) );
        BERvector = transpose ( statsMatrix(: , : ,2) ./ statsMatrix(: , : ,3) );        
        BERvector = BERvector ./ berFix;
        BERvector = BERvector(: , 2:maxN );
        
        % Plot SEP
        figure
        semilogy(snrvalues,SEPvector(: , 2:maxN ),'-*')        
        legend("n="+string(2:maxN)+" "+"",'Location','Best')
        xlabel('SNR (dB)') 
        ylabel('Symbol Error Probability') 
        title("SEP versus received SNR")
        grid on

        % Plot BER
        figure
        semilogy(snrvalues,BERvector,'-v') 
        legend("n="+string(2:maxN)+" "+"",'Location','Best')        
        xlabel('SNR (dB)') 
        ylabel('Bit Error Probability') 
        title("BER versus received SNR")
        grid on
else
        fprintf('\nThere were no files in the current directory to load data\n');
        
end

fprintf('\n');
        
