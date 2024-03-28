
% This runs the simulations to calculate the SEP and BER 
maxN = 10 ;
maxSNR = 40 ;
totalSymbolsSent = 30 ; % This is the power of 2 of the symbols to sent  2^totalSymbolsSent
                        % If it left 30 for low SNR values where the SEP is very low it will take a lot of time !!!
saveResults = true; 
%saveResults = false; 

% This is the Date that it was run
dateAndTimeNow = datetime('now');
DateNumber = datenum(dateAndTimeNow) ; % Convert to double
datetime(DateNumber,'ConvertFrom','datenum') ; % Convert it to datetime from double

thisRunStatsMatrix = zeros ( maxN, maxSNR, 1 , 7 );
%  n , snr , Number of Run ,  vector of data
% vector of data :  = [ SEP BER NOS Time DetecAlgoID Who When ] 


% The name of the file to save the results 
%nameOfThePointsFile = 'MultiRunsValuesForSEP&BER_Nick.mat' ;
%nameOfThePointsFile = "MultiRunsValuesForSEP&BER_Katerina.mat" ;
nameOfThePointsFile = 'MultiRunsValuesForSEP&BER_Extra.mat' ;

% Also the 6th element of the array around the line 110 should be changed to the appropriate integer 


distance = 2;
snr = 1:1:maxSNR; 
SERvector = zeros(length(snr),maxN);
BERvector = zeros(length(snr),maxN);


tStart = tic; % Pair 1 : tic 

for currentN = 2:maxN
        %currentN
        
        m = 2^currentN;                         %number of symbols in the constellation
        t = 2^(totalSymbolsSent-currentN); 	%number of symbols sent
        
        if rem(currentN,2)==0 && currentN>=2
            p = 1;
            a = 2^(currentN/2); 
        elseif rem(currentN,2)==1 && currentN>3
            p = 2;
            a = 2^((currentN-3)/2);
        elseif currentN==3
            p = 3;
            a = 0; % Some random value to not crash when the DetectionFast is called, it is not usefull since n=3 is a special case
        end
        
        
        [SymbolCoordinates,SymbolData,constellationVector , constellationGrayCodeVector ,constPower] = RegularHQAM(currentN,distance);
        sigpower = pow2db(constPower);
        
        IndexOfSNR =1 ;
        NumberOfSymbolsSent= t*m; % = 2 ^ totalSymbolsSent ; 
        
        for currentSNR = snr
                
		%{     
                if statsMatrix(currentN , currentSNR , 1) > 1000 || statsMatrix(currentN , currentSNR , 1)==0
                        % This requiaries to have already run the Final_SEP_BER_graph_of_all_data.m to create the statsMatrix 
                        % It skips the points that have already a lot of symbol errors which means the graph should be smooth there
                        % It is used to smooth out mostly the last points of the diagrams where the SEP is very low and the expected errors are very little
                        fprintf('\n Skipped  N= %d,\tSNR = %d ' , currentN , currentSNR );
                        continue
                end
                %}
                
                BER=0;
                SER=0;                
                
                %currentSNR                
                
                tStartForNOS = tic; % Pair 2 : tic 
                        
                
                parfor k=1:t 
                    
                    % We create a vector with all the symbols of the constellation with noise 
                    symbolsWithNoise = awgn(constellationVector,currentSNR,sigpower);
                    for i=1:m
                            % Depending on the detection algorithm the appropriate change must be done to the storing of the data in the array at line 110
                            %[SERtemp , BERtemp ] =  DetectionLinearSlow(constellationVector , constellationGrayCodeVector , symbolsWithNoise(i) , constellationGrayCodeVector(i)  , m , distance  );
                            [SERtemp , BERtemp ] = DetectionFast(SymbolCoordinates , SymbolData , symbolsWithNoise(i) , constellationGrayCodeVector(i)  , p , distance , a);
                            
                            SER = SER+ SERtemp;
                            BER = BER+ BERtemp;                            

                    end
                    
                end 
               

                % Calculate the BER and SER for the currentN and currentSNR
                BERper=BER/NumberOfSymbolsSent;
                SERper=SER/NumberOfSymbolsSent;
                
                % Store them to plot them 
                SERvector(IndexOfSNR, currentN) = SERper;
                BERvector(IndexOfSNR, currentN) = BERper;
                IndexOfSNR = IndexOfSNR +1 ;
                
                % The time it took to calculate the BER and SER for the currentN and currentSNR
                tEndForNOS = toc(tStartForNOS); % Pair 2 : toc
                
                % Store all the data to combine them from different runs 
                thisRunStatsMatrix( currentN, currentSNR, 1 , : )= [SER BER NumberOfSymbolsSent tEndForNOS 20 3 DateNumber];
                %  n , snr , Number of Run ,  vector of data
                % vector of data :  = [ SEP BER NOS Time DetecAlgoID Who When ] 
                
                %{
                DetecAlgoID :  
                        10 DetectionLinearSlow
                        20 DetectionFast
                Who:
                        1 Nick          odd     Laptop
                        2 Katerina      even    PC
                        3 Extra         ???     ???
                        4 ???           ???	???
                                                
                %}
                
                fprintf('\n Completion: %f%% , Time for this point: %f sec,\t Total time: %f sec,\t\t N= %d,\tSNR = %d ' , (9* (currentN-2)/(maxN-2) +  currentSNR/maxSNR ) *10 , tEndForNOS , toc(tStart) ,currentN , currentSNR );
                
                
                if ( SERper <  1 /NumberOfSymbolsSent )
                        %SERper
                        % It can stop one loop earlier (before we have a zero)
                        % This stops the loop for high SNR when the expected BER would be too low to evaluate whith the NumberOfSymbolsSent 
                        break
                end
                
                
                
        end        
end

tEnd = toc(tStart); % Pair 1 : toc               % Total time program was running 

% Plot SER 
figure
semilogy(snr,SERvector,'-*')


% Plot BER
figure
semilogy(snr,BERvector,'-v')



fprintf('\n\n Succesfull completion! \n Total time: %d min and  %f sec \n' , floor(tEnd/60) ,  mod(tEnd , 60) );
sound(sin(1:3000));     % Beep when it's over       

if  exist( nameOfThePointsFile , 'file') == 2 && saveResults
        %This means the file exists :) 
        load(nameOfThePointsFile);
        fprintf('\nLoaded the file named : %s with last SEP and BER values \n' , nameOfThePointsFile);
        %nameOfTheFile % just to print the name of the file 
        
        % Now we need to add the new data points to the old one 
        i = 1 + size(MultiRunsValues, 3);
        MultiRunsValues(:, :, i , : )  = thisRunStatsMatrix(:, :, 1  , : );
       
        save(nameOfThePointsFile,'MultiRunsValues');
        
elseif saveResults
       
        MultiRunsValues = thisRunStatsMatrix ;
        save(nameOfThePointsFile,'MultiRunsValues'); 
        
end












