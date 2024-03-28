


Final_SEP_BER_graph_of_all_data;
SEPsimulations = SEPvector(: , 2:10);  % Delete the values for n=1 to plot them right
clearvars -except SEPsimulations BERvector statsMatrix statsMatrixLinearSlow statsMatrixFastDetection    %deletes all variables except X in workspace


%diagramPrints = [ false false false false false false false false false false false false false false];
diagramPrints = [ true true  true true  true true  true true  true true  true true  true true  true true ];

printIndex=0; 


b = [4 8 10 13 22 27 46 64 128 ] ; % Thrassos values
%b = [4 8 12 18 28 42 60 90 124] ;  % Our values

m = [ 4   8   16   32   64   128  256  512  1024 ]; % Total number of symbols in constellation 
%b % Total number of external symbols in the constellation
%r % Special XD
%k 

SNRvalues = 1:40 ; % The snr we want to calculate the SEP 

Es = [ 
        2.0000
        4.5000
        9       % n = 4
        17.7500
        37
        72
        149
        289
        597
        %1157
        %2389
        %4629
        %9557
        %18517
        %38229
        %74069
        %152917
        %296277
].';

k = [   0.9999999
        0.9711505
        0.8711505       % n = 4
        0.7233274
        0.5222431
        0.5088351
        0.3936315
        0.3672311
        0.2982858
        ].';

kupperBound = zeros( 1 , 9); 




dmin = 2 ;      % Anything we want
R = dmin / 2;
Rtonos = dmin / sqrt(3);

r = k * Rtonos + (1-k) * R ;
rupperBound = kupperBound * Rtonos + (1-kupperBound) * R ;

%SNR =  10 log Es/ No ; 
%No = Es / 10^(snr / 10 )  ;       % This determines the SNR (with dmin) 
%SEPapr = (2*m-b) ./ (2*m) .* exp( - r.^2 ./ No) + b./m .* qfunc( r.* sqrt(2./No));


SEPaprox = zeros ( 40 , 9) ; 
SEPaproxUpperBound = zeros ( 40 , 9) ;

for snr = SNRvalues
        
        No = Es / 10^(snr / 10 )  ;       % This determines the SNR (with dmin) 
        %r = k * Rtonos + (1-k) * R ; 
        
        SEPapr = (2*m-b) ./ (2*m) .* exp( - r.^2 ./ No) + b./m .* qfunc( r.* sqrt(2./No));        
        SEPaprUpperBound = (2*m-b) ./ (2*m) .* exp( - rupperBound.^2 ./ No) + b./m .* qfunc( rupperBound.* sqrt(2./No));
        
        SEPaprox( snr , :) = SEPapr ;
        SEPaproxUpperBound( snr , :) = SEPaprUpperBound ;
end 

clear SEPapr SEPaprUpperBound ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



snrvalues = 1:40;
SEPvector = SEPsimulations;

% BER aproximation using SEP and Gp

% The name of the file we save the values 
nameOfTheGpFile = 'GpValues.mat' ;
if  exist( nameOfTheGpFile , 'file') == 2  
        %This means the file exists :) 
        load(nameOfTheGpFile);
        fprintf('\nLoaded the file named : %s with pre-calculated values for the Gp \n' , nameOfTheGpFile);
        %nameOfTheFile % just to print the name of the file 
else 
        GpValues = [0;1.16666666666667;1.27500000000000;1.23750000000000;1.38854166666667;1.28229166666667;1.36354166666667;1.30703125000000;1.35156250000000;1.31998697916667;1.34322916666667;1.32661132812500;1.33847656250000;1.32996012369792;1.33595377604167;1.33164367675781;1.33465576171875;1.33248774210612;1.33399759928385;1.33291034698486];
end

berFixDueToN = zeros(40,9);
berFixDueToGp = zeros(40,9);

for i = snrvalues
        berFixDueToN(i,:) = 2:10;
        berFixDueToGp(i,:) = GpValues(2:10,1);
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot Relative Error of BER Approximations (from BER = SEP*Gp/n) from Simulations

BERaprox = SEPvector .* berFixDueToGp  ./ berFixDueToN;
RelativeBERError = abs(BERvector - BERaprox ) ./ BERvector;

legendRelativeBER_Error= [ "n="+string(2:9)+"  Gp="+string(GpValues(2:9,1).') "n="+string(10:10)+" Gp="+string(GpValues(10:10,1).')];
printIndex = printIndex + 1; 
if diagramPrints(printIndex)          
        
        figure
        semilogy(snrvalues,RelativeBERError,'-*') 

        
        legend(legendRelativeBER_Error,'Location','Best') ;
        xlabel('SNR (dB)') 
        ylabel('Relative BER Error') 
        title("Relative BER Error of Approximations from Simulations ")
        grid on
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot SEP with Approximations and Upper Bounds

legendONE =  "n="+string(2:10)+" " ;

printIndex = printIndex + 1; 
if diagramPrints(printIndex)        
       
        figure
        semilogy(snrvalues,SEPsimulations,'-*') 
        hold on 
        semilogy(snrvalues,SEPaprox,'-.o')  
        hold on 
        semilogy(snrvalues, SEPaproxUpperBound,'--s') 

        ylim([10^(-9) 1])
        
        legendAll = [ legendONE  legendONE legendONE] ;

        lg = legend(legendAll,'Location','Best','NumColumns',3) ;
        title(lg, "Simulations" + "   " + "Approximation" + "   " + "Upper Bound" )

        xlabel('SNR (dB)') 
        ylabel('Symbol Error Probability') 
        title("SEP versus received SNR for Simulations, Approximations and Upper Bounds")
        grid on
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Plot BER with Approximations 
printIndex = printIndex + 1; 
if diagramPrints(printIndex)        
        
        figure
        semilogy(snrvalues,BERvector,'-*') 
        hold on 
        semilogy(snrvalues,BERaprox,'-.o')  
        legendBERgraph = [legendONE legendRelativeBER_Error ] ;
        %ylim([10^(-9) 1])
        lgBER = legend(legendBERgraph,'Location','Best','NumColumns',2) ;
        title(lgBER, "Simulations" + "           " + "Approximation   "  )

        xlabel('SNR (dB)') 
        ylabel('BER') 
        title("BER versus received SNR for Simulations and Approximations")
        grid on
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot Relative Error of SEP Approximations from Simulations

RelativeSEPError = abs(SEPsimulations - SEPaprox ) ./ SEPsimulations;

printIndex = printIndex + 1; 
if diagramPrints(printIndex)        
        
        legendRelativeError=  "n="+string(2:10)+" k="+string(k) ;
        figure
        semilogy(snrvalues,RelativeSEPError,'-*') 

        legend(legendRelativeError,'Location','Best') ;
        xlabel('SNR (dB)') 
        ylabel('Relative Error') 
        title("Relative Error of SEP Approximations from Simulations")
        grid on
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Plot the difference of the Upper Bound for SEP from the Simulations (It should be positive if not Matlab throws a warning)

UpperBoundCheck = SEPaproxUpperBound- SEPsimulations ;

printIndex = printIndex + 1; 
if diagramPrints(printIndex)        
       
        figure
        semilogy(snrvalues,UpperBoundCheck,'-*') 

        legend(legendONE,'Location','Best') ;

        ylim([10^(-9) 1])

        xlabel('SNR (dB)') 
        ylabel('Absolute Error') 
        title("Absolute Error of Upper Bound from Simulations")
        grid on
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot the relative error of the Upper Bound for SEP from the Simulations 

UpperBoundRelativeError = UpperBoundCheck ./ SEPsimulations ; 

printIndex = printIndex + 1; 
if diagramPrints(printIndex)        
       
        figure
        semilogy(snrvalues,UpperBoundRelativeError,'-*') 

        legend(legendONE,'Location','Best') ;

        ylim([10^(-4) 10])

        xlabel('SNR (dB)') 
        ylabel('Relative Error') 
        title("Relative Error of Upper Bound from Simulations")
        grid on
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Averege time to detect per Symbol 

%statsMatrixLinearSlow ;
%statsMatrixFastDetection;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Average Time to detect a Symbol with FastDetection

AverageTimePerSymbolFastDetection = statsMatrixFastDetection(:,:, 4 ) ./ statsMatrixFastDetection(:,:, 3 ) ;
AverageTimePerSymbolFastDetection = AverageTimePerSymbolFastDetection( 2:10,: );  % Delete the values for n=1 to plot them right

% vector of data :  = [ SEP BER NOS Time DetecAlgoID Who When ] 
printIndex = printIndex + 1; 
if diagramPrints(printIndex)        
        
        figure
        semilogy(snrvalues,AverageTimePerSymbolFastDetection,'-*') 

        legend(legendONE,'Location','Best') ;

        ylim([10^(-9) 1])

        xlabel('SNR (dB)') 
        ylabel('Average Time per Symbol (sec)') 
        title("Average Time to detect a Symbol with FastDetection")
        grid on
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Average Time to detect a Symbol with LinearSlow

AverageTimePerSymbolLinearSlow = statsMatrixLinearSlow(:,:, 4 ) ./ statsMatrixLinearSlow(:,:, 3 ) ;
AverageTimePerSymbolLinearSlow = AverageTimePerSymbolLinearSlow( 2:10,: );      % Delete the values for n=1 to plot them right

printIndex = printIndex + 1; 
if diagramPrints(printIndex)        
       
        figure
        semilogy(snrvalues,AverageTimePerSymbolLinearSlow,'-*') 

        legend(legendONE,'Location','Best') ;

        ylim([10^(-9) 1])

        xlabel('SNR (dB)') 
        ylabel('Average Time per Symbol (sec)') 
        title("Average Time to detect a Symbol with LinearSlow")
        grid on
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Average Time to detect a Symbol for ALL the Simulation data

%AverageTimePerSymbolTotal = statsMatrix(:,:, 4 ) ./ statsMatrix(:,:, 3 ) ;
%AverageTimePerSymbolTotal = AverageTimePerSymbolTotal( 2:10,: );  % Delete the values for n=1 to plot them right




