%% Spectral characteristics of data

%% Load
load ./data/alignedEMG.mat

%% Select dataset:
subject=2;
epoch=1;
%% Decompose average and trial-to-trial
data=alignedEMG{subject,epoch}.Data;
avg=mean(data,3);
dc=mean(avg,1);
T2T=permute(data-avg,[2,1,3]);
T2T=T2T(:,:)';

T2T=T2T./sqrt(sum(T2T.^2));
fullData=permute(data-dc,[2,1,3]);
fullData=fullData(:,:)';
fullData=fullData./sqrt(sum(fullData.^2));

%% Spectral characteristics visualization
figure
N=size(fullData,1);
subplot(8,4,1:4:10) %Autocorrelation matrix full process
C=fullData*fullData';
imagesc(C); caxis(norm(fullData,'fro')^2*[-1 1]/N)
title(['Centered, PR=' num2str(PReff(C))])
subplot(8,4,13)
plot(C(1,:))
subplot(8,4,[2:4:10])
Ct=cyclicToeplitize(C);
imagesc(Ct); caxis(norm(fullData,'fro')^2*[-1 1]/N)
title(['Toeplitz, PR=' num2str(PReff(Ct))])
subplot(8,4,14)
plot(Ct(1,:))
subplot(8,4,3:4:12) %Autocorrelation T2T
C=T2T*T2T';
imagesc(C); caxis(norm(T2T,'fro')^2*[-1 1]/N)
title(['T2T, PR=' num2str(PReff(C))])
subplot(8,4,4:4:12)
Ct=toeplitize(C);
imagesc(Ct); caxis(norm(T2T,'fro')^2*[-1 1]/N)
title(['Toeplitz, PR=' num2str(PReff(Ct))])
subplot(8,4,16)
plot(Ct(1,:))

subplot(2,2,3) %FFT avg & T2T
plot(1e5*mean(abs(fft(fullData)),2).^2);
hold on
plot(50e5*mean(abs(fft(T2T)),2).^2);
axis([1 300 0 1])

subplot(2,2,4) %loglog FFT
l1=loglog(1e5*mean(abs(fft(fullData)),2).^2,'o','MarkerSize',8,'DisplayName','Periodic component');
l1.MarkerFaceColor=l1.Color;
hold on
l2=loglog(50e5*mean(abs(fft(T2T)),2).^2,'DisplayName','Trial-to-trial variability');
axis([10 500 2e-2 2])
grid on
for i=2%1:3
loglog([30,300],1*[1 .1]/i^2,'Color',[255,105,180]/255,'LineWidth',2)
end
text(15,1/2,'f^{-1}','Color',[255,105,180]/255,'Fontsize',10,'FontWeight','bold')

for i=1%:3
loglog([30,300],1*[1 .0316]/i^2,'Color',[160,82,45]/255,'LineWidth',2)
end
text(15,(1.5),'f^{-1.5}','Color',[160,82,45]/255,'Fontsize',10,'FontWeight','bold')
legend([l1 l2])