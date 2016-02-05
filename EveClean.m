%% Even-Skipped Analysis
%% Load Data
load 2x;
Data = Data([1:11])
load Background;
Factor = 0.6431;
EL=[Data.EL];  
[sEL,seq]=sort(EL);
aEL=EL*Factor;
disp('the sample number is');length(aEL)
disp('the mean EL is');mean(aEL)
disp('the std of EL is');std(aEL)
Bins = length(Data(1).relx);
L = floor(EL/Bins);
xx = 1/Bins/2:1/Bins:1-1/Bins/2;
L2 = round(900/Bins);
C = rand(length(Data),3);
%%
boundary=zeros(length(Data),24);
figure, 
for i=1:length(Data), 
    x=Data(i).relx;
    y=Data(i).relHb;
    x=x/EL(i);
    y2=sort(y);                           
    thresholdmax=mean(y2(length(y2)-4:length(y2)));
    thresholdmin=mean(y2(1:10));
%     yy=(y-thresholdmin)/(thresholdmax-thresholdmin); %Switch Normalization
    yy=(y-thresholdmin) %Normalization off
    scatter(x,yy,[],C(i,:),'.');hold on
    PosteriorHb(i)=max(yy(find(x>0.65)));
    AnteriorHb(i)=max(yy(find(x<0.32)));
    yMinusHalf=yy(:)-0.5;
    yShift1=[0; yMinusHalf];yShift2=[yMinusHalf;0];
    yShiftMinus=yShift1.*yShift2;
    foundHalf=find(yShiftMinus<0);
    for q=1:length(foundHalf)
        if yy(foundHalf(q)-1)>yy(foundHalf(q))
            boundary(i,q)=interp1(yy(foundHalf(q)-1:foundHalf(q)),x(foundHalf(q)-1:foundHalf(q)),0.5);
        end
    end
    yyy(i,:)=yy;
    sl(i,:)=gradient(yy);
    ssl(i,:)=gradient(smooth(yy));
end
