%% Load data
load('../data/alignedEMG.mat') %alignedEMG is the variable name
load('../data/EMG1.mat') %Non-aligned BASELINE data
%% Get example data:
subj=5;
%Grab aligned data, and transform to 2D array: time x muscles
alignedData=permute(alignedEMG{subj,1}.Data,[1,3,2]); 
alignedData=reshape(alignedData,size(alignedData,1)*size(alignedData,2),size(alignedData,3));

%Grab the SAME time interval, unaligned
t0=alignedEMG{subj,3}.eventTimes(1,1);
t1=alignedEMG{subj,3}.eventTimes(1,end);
unalignedData=EMG3{subj}.Data(EMG3{subj}.Time>=t0 & EMG3{subj}.Time<=t1,:); %This is 1000 samples longer than the aligned version for subj=5, meaning there is slight downsampling when aligning, shouldn't matter since the signal is sampled at 1kHz

%% Pre-process:
%Column normalization:
uD=unalignedData;
aD=alignedData;
uD=unalignedData./sqrt(sum(unalignedData.^2,2));
aD=alignedData./sqrt(sum(alignedData.^2));
M=size(uD,2);
%% Spectrum estimation:
Fu=fft(uD,[],1);
Fa=fft(aD,[],1);
su=sum(abs(Fu).^2,2);
sa=sum(abs(Fa).^2,2);
Ca=dct(aD,[],1);
ta=sum(abs(Ca).^2,2);
Da=fcst(aD);
wa=sum(abs(Da).^2,2);
%% vaps!
lu=eig(uD'*uD);
la=eig(aD'*aD);

%% Figure
figure;
subplot(1,2,1)
hold on
%plot(VAF(lu),'DisplayName','UnalignedVAF')
plot(VAF(la),'DisplayName','AlignedVAF')
%plot(VAF(su),'DisplayName','Unaligned PSD bound')
c=VAF(sa);
plot(c,'DisplayName','Aligned PSD bound')
plot(c(1:M)+(1-c(1:M)).*[1:M]'/M,'DisplayName','Experimental Aligned PSD bound')
d=VAF(ta);
%plot(d,'DisplayName','Aligned DCT bound')
%plot(d(1:M)+(1-d(1:M)).*[1:M]'/M,'DisplayName','Experimental Aligned DCT bound')
e=VAF(wa);
plot(e,'DisplayName','Aligned FCST bound')
plot(e(1:M)+(1-e(1:M)).*[1:M]'/M,'DisplayName','Experimental Aligned FCST bound')
axis([1 M 0 1])
legend('Location','SouthEast')

subplot(1,2,2)
hold on
title('Aligned VAF excess')
%plot(VAF(lu),'DisplayName','UnalignedVAF')
%plot(VAF(la),'DisplayName','AlignedVAF')
%plot(VAF(su),'DisplayName','Unaligned PSD bound')
c=VAF(sa);
plot(VAF(la)- c(1:M),'DisplayName','PSD bound excess')
plot(VAF(la)-c(1:M)-(1-c(1:M)).*[1:M]'/M ,'DisplayName','Experimental PSD bound excess')
%plot(d,'DisplayName','Aligned DCT bound')
%plot(d(1:M)+(1-d(1:M)).*[1:M]'/M,'DisplayName','Experimental Aligned DCT bound')
plot(VAF(la)-e(1:M),'DisplayName','FCST bound excess')
plot(VAF(la)-e(1:M)-(1-e(1:M)).*[1:M]'/M,'DisplayName','Experimental Aligned FCST excess')
axis([1 M 0 .2])
legend('Location','NorthEast')
