%HQAM 

%inputs
%distance = sqrt(3);   	%minimum distance between 2 symbols
distance = 2;
n = 5;                	%order of the constellation
m= 2^n;
c = 2; % do we need this ? 

[SymbolCoordinates,SymbolCoordinates2,SymbolCoordinates2Transpose,SymbolData,refConst,constPower] = RegularHQAM(n,distance);


t = 100;                %number of symbols sent
sigpower = pow2db(mean(abs(refConst).^2));
%sigpower = constPower;
%rng('shuffle');
h = scatterplot(refConst,[],[],'r*');
grid
drawnow
hold on
for i=1:t
    x = randi(m);
    rxSig = awgn(refConst(x),15,sigpower);
    scatterplot(rxSig,[],[], 'b.',h)
    drawnow
end 

hold off




%Construct a constellation diagram
constDiagram = comm.ConstellationDiagram('Title',string(m)+ ' HQAM', ...
    'XLimits',[-sqrt(m)*distance/2-c*distance sqrt(m)*distance/2+c*distance],'YLimits',[-sqrt(m)*distance/2 sqrt(m)*distance/2], ...
    'ReferenceConstellation',refConst, ...
    'ReferenceMarker','*','ReferenceColor',[0 1 0]);

%Plot the customized constellation
constDiagram(refConst)