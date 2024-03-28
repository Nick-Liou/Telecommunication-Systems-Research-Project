
%Testing the new Fast detection 



maxN = 8;
distance = 2;
step = 0.05 ;  % This is how fine the grid of points that tests the DetectionFast is 

edgeCase = zeros(1,maxN);
wrongDetection = zeros(1,maxN);

tStart = tic; % Pair 1 : tic 


for n = 2:maxN        
        
        tStartForN = tic; % Pair 2 : tic 
        
        m= 2^n;
                
        if rem(n,2)==0 && n>=2
            p = 1;
            a = 2^(n/2);         
            testLength = a*distance /2 ;    
        elseif rem(n,2)==1 && n>3
            p = 2;
            a = 2^((n-3)/2);
            testLength = 3*a*distance /2 ;  
        elseif n==3
            p = 3;
            a = 2^((n-3)/2);
            testLength = 3*a*distance /2 ;              
        end
        
                        
        [SymbolCoordinates,SymbolData,constellationVector, constellationGrayCodeVector,~] = RegularHQAM(n,distance);
        
        % This makes the testing region smaller to run faster 
        %testLength = testLength *1/3;
        
        
        squareDots = -testLength: step :testLength;

        
        h = scatterplot(constellationVector,[],[],'r*');
        grid
        drawnow
        hold on
        
        for x = 1:length(squareDots)
                for y = 1:length(squareDots)
                        point = squareDots(x) + 1i*squareDots(y);
                        
                        %Now we have the point and we need to test
                        
                        
                        
                        % This needs a different version of the function to work (which has extra logic and computations and thus is slower)
                        [~,~, slowGrayCode] = DetectionLinearSlowTesting(constellationVector , constellationGrayCodeVector , point  , -2   , m , distance  ); 
                        % it must be changed a bit to return the GrayCode of the symbol
                        
                        [~ , ~ , fastGrayCode ] = DetectionFastTesting(SymbolCoordinates , SymbolData , point , slowGrayCode  , p , distance , a);
                        
                        if fastGrayCode == -2 
                                edgeCase(n) = edgeCase(n) +1 ; 
                                % This prints a green dot if it is one of the outside symbols that the old DetectionFast could not handle 
                                % Now since for those symbols it calls the DetectionLinearSlow2D it gets correctly 
                                % To see what symbols the DetectionFastTesting does not find we need to add in the final lines
                                % of the DetectionFastTesting function where it calls the DetectionLinearSlow2D the GrayCode = -2
                                scatterplot(point,[],[], 'g.',h)
                                
                        elseif slowGrayCode ~= fastGrayCode 
                                wrongDetection(n) = wrongDetection(n) +1 ;
                                if fastGrayCode == -1
                                        % This print a magenta dot if the coordinates of the symbol detected from the DetectionFastTesting correspond to a symbol
                                        % that could be stored in the 2D matrix that the constellation symbols are stored but in that cell it does not have a symbol
                                        % This does not happen 
                                        scatterplot(point,[],[], 'm.',h)
                                else
                                        % This prints a blue dot if the symbol detected from the DetectionLinearSlowTesting and DetectionFastTesting
                                        % do not match, which can happen if the DetectionFastTesting detects a wrong symbol (which it does not) 
                                        % and if the symbol is on the line of the decision boundaries (which happens ! ) 
                                        scatterplot(point,[],[], 'b.',h)
                                        
                                end
                        else 
                                % This would print a red dot if the DetectionLinearSlowTesting and DetectionFastTesting find the same symbol
                                % This can be used to see how fine the grid is and otherwise is commented 
                                % for visual clarity and speed
                                %scatterplot(point,[],[], 'r.',h)
                        end
                        
                        %drawnow  % This would print the dots one by one to be able to see the progress but it would make it much much slower 
                        

                end
                %drawnow  % This would print the dots line by line to be able to see the progress but it would be slower

        end
        drawnow
        
        % How much time it took to run for each n         
        tEndForN = toc(tStartForN); % Pair 2 : toc        
        fprintf('\n Total time to test for n = %d was %f sec' , n , tEndForN );

        
end


hold off

tEnd = toc(tStart); % Pair 1 : toc               % Total time program was running 


fprintf('\n Total time to test for all n with step %d was : %d min and  %f sec \n' , step ,  floor(tEnd/60) ,  mod(tEnd , 60) );













