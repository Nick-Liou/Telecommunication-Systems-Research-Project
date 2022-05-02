


n = 10;
distance = 2; 
%snr = 35 ; 
snr = 0:3:20; 
SERvector = zeros(length(snr),n);
BERvector = zeros(length(snr),n);

tic
for currentN=2:n

        m = 2^currentN;                 %number of symbols
        t = 2^(15-currentN);           	%number of symbols sent
        
        [~,~,~,~,constellationVector, constellationGrayCodeVector,constPower] = RegularHQAM(currentN,distance);
        sigpower = pow2db(constPower);
        
        IndexOfSNR =1 ;
        NumberOfSymbolsSent= t*m;
        %tic
        for currentSNR = snr
                %currentSNR
                
                BER=0;
                SER=0;
                

                
                parfor k=1:t
                    %tic
                    %x = randi(m);
                    rxSig = awgn(constellationVector,currentSNR,sigpower);
                    for i=1:m
                            [SERtemp , BERtemp ] =  DetectionLinearSlow(constellationVector , constellationGrayCodeVector , rxSig(i) , constellationGrayCodeVector(i)  , m , distance  );
                            SER = SER+ SERtemp;
                            BER = BER+ BERtemp;
                            

                    end
                    %NumberOfSymbolsSent= NumberOfSymbolsSent + m;


                end 
               


                BERper=BER/NumberOfSymbolsSent;
                SERper=SER/NumberOfSymbolsSent;
                
                SERvector(IndexOfSNR, currentN) = SERper;
                BERvector(IndexOfSNR, currentN) = BERper;
                IndexOfSNR = IndexOfSNR +1 ;
                
                if ( SERper <  1 /NumberOfSymbolsSent )
                        %SERper
                        % It can stop one loop earlier (before we have a zero) Tho now its for 3 dB 
                        % This stops the loop for high SNR when the expected BER would be too low to evaluate whith the NumberOfSymbolsSent 
                        break
                end
        end
        toc
end


figure
%semilogy(snr,SERvector(:,2),'g',snr,SERvector(:,3),'b--o',snr,SERvector(:,4),'c*')

semilogy(snr,SERvector,'-*')



%{
[~,~,~,~,constellationVector, constellationGrayCodeVector,constPower] = RegularHQAM(n,distance);

sigpower = pow2db(constPower);

rng('shuffle');
h = scatterplot(constellationVector,[],[],'r*');
grid
drawnow
hold on


% DetectionLinearSlow(constellationVector , constellationGrayCodeVector , symbolsWithNoise , grayCodeOfSymbolSent  , m , d  ) 

BER=0;
SER=0;
NumberOfSymbolsSent=0;

if n~=3 || true
        tic
        parfor k=1:t
            %tic
            x = randi(m);
            rxSig = awgn(constellationVector,snr,sigpower);
            for i=1:m
                    [SERtemp , BERtemp ] =  DetectionLinearSlow(constellationVector , constellationGrayCodeVector , rxSig(i) , constellationGrayCodeVector(i)  , m , distance  );
                    SER = SER+ SERtemp;
                    BER = BER+ BERtemp;
                        
            end
            NumberOfSymbolsSent= NumberOfSymbolsSent + m;
            
            
            %scatterplot(rxSig,[],[], 'b.',h)
            
            %drawnow
            %toc
        end 
        toc
        
end
%{
else
         for k=1:t
            %tic
            x = randi(m);
            rxSig = awgn(complex(constellationVector(x)),snr,sigpower);
            scatterplot(rxSig,[],[], 'b.',h)
            drawnow
            %toc
        end 
end
%}





drawnow
hold off


BERper=BER/NumberOfSymbolsSent;
SERper=SER/NumberOfSymbolsSent;

%}














