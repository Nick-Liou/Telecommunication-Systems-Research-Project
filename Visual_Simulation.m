%HQAM 

%minimum distance between 2 symbols
distance = 2;
n = 4;                	%power of the constellation
m= 2^n;

SNR = 15 ;              % SNR
NOS = 30;                %number of symbols sent


[~,~,constellationVector, ~,constPower] = RegularHQAM(n,distance);

%sigpower = pow2db(mean(abs(constellationVector).^2));
sigpower = pow2db(constPower);

h = scatterplot(constellationVector,[],[],'r*');
grid
drawnow
hold on
for i=1:NOS
    x = randi(m);
    rxSig = awgn(constellationVector(x),SNR,sigpower);
    scatterplot(rxSig,[],[], 'b.',h)
    drawnow % draws each symbol sent one by one 
end 

drawnow
hold off




c = 2;
%Construct a constellation diagram
constDiagram = comm.ConstellationDiagram('Title',string(m)+ ' HQAM', ...
    'XLimits',[-sqrt(m)*distance/2-c*distance sqrt(m)*distance/2+c*distance],'YLimits',[-sqrt(m)*distance/2 sqrt(m)*distance/2], ...
    'ReferenceConstellation',constellationVector, ...
    'ReferenceMarker','*','ReferenceColor',[1 0 0]);

%Plot the customized constellation
constDiagram(constellationVector)