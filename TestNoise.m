


n = 4;
distance = 2; 
snr = 10 ;
m = 2^n;                 %number of symbols
t = 50;                %number of symbols sent

[~,~,~,~,refConst,constPower] = RegularHQAM(n,distance);

sigpower = pow2db(constPower);

rng('shuffle');
h = scatterplot(refConst,[],[],'r*');
grid
drawnow
hold on
for i=1:t
    %tic
    x = randi(m);
    rxSig = awgn(refConst(x),snr,sigpower);
    scatterplot(rxSig,[],[], 'b.',h)
    %drawnow
    %toc
end 

drawnow
hold off