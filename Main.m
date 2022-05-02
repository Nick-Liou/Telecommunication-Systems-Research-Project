%HQAM 

%inputs
%distance = sqrt(3);   	%minimum distance between 2 symbols
distance = 2;
n = 4;                	%order of the constellation
m= 2^n;
c = 2; % do we need this ? 

[SymbolCoordinates,SymbolCoordinates2,SymbolCoordinates2Transpose,SymbolData,constellationVector, constellationGrayCodeVector,constPower] = RegularHQAM(n,distance);


t = 10;                %number of symbols sent
sigpower = pow2db(mean(abs(constellationVector).^2));
%sigpower = constPower;
%rng('shuffle');
h = scatterplot(constellationVector,[],[],'r*');
grid
drawnow
hold on
for i=1:t
    x = randi(m);
    rxSig = awgn(constellationVector(x),15,sigpower);
    scatterplot(rxSig,[],[], 'b.',h)
    drawnow
end 

hold off




%Construct a constellation diagram
constDiagram = comm.ConstellationDiagram('Title',string(m)+ ' HQAM', ...
    'XLimits',[-sqrt(m)*distance/2-c*distance sqrt(m)*distance/2+c*distance],'YLimits',[-sqrt(m)*distance/2 sqrt(m)*distance/2], ...
    'ReferenceConstellation',constellationVector, ...
    'ReferenceMarker','*','ReferenceColor',[0 1 0]);

%Plot the customized constellation
constDiagram(constellationVector)