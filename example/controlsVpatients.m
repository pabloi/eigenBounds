%Compare stroke and controls
load('../data/alignedEMGControls.mat') %alignedEMG is the variable name
load('../data/alignedEMGPatients.mat') %alignedEMG is the variable name

%%
l=nan(32,30);
v=nan(32,30);
s=nan(32,7296);
t=nan(32,7296);
for i=1:16
    for k=1:2
        switch k
            case 1
                D=alignedEMGControls{i,1}.Data; %Baseline walking
            case 2
                D=alignedEMGPatients{i,1}.Data; %Baseline walking
        end
            D=D(:,:,1:19);
            %D=D-mean(D,3); %Mean across strides
            D=permute(D,[1,3,2]);
            D=reshape(D,size(D,1)*size(D,2),size(D,3));
            %D=D-mean(D);
            D=D./sqrt(sum(D.^2,1));
            l(i+16*(k-1),:)=eig(D'*D);
            v(i+16*(k-1),:)=VAF(l(i+16*(k-1),:));
            
            F=fft(D,[],1);
            aux=VAF(sum(abs(F).^2,2));
            s(i+16*(k-1),:)=aux;
            F=fcst(D);
            aux=VAF(sum(abs(F).^2,2));
            t(i+16*(k-1),:)=aux;
    end
end

%% 
v1=v(1:16,:);
%v2=v(16+[1:5,7:10,12:16],:);
v2=v(17:32,:);
s1=t(1:16,1:30);
%s2=s(16+[1:5,7:10,12:16],1:30);
s2=t(17:32,1:30);
for i=1:size(v1,1)
   s1(i,:)=improveBound(s1(i,:),30); 
   s2(i,:)=improveBound(s2(i,:),30); 
end
for i=1:30
    [hh(i),pp(i)]=ttest2(v1(:,i),v2(:,i));
    %[pp(i),hh(i)]=ranksum(v1(:,i),v2(:,i));
    [hh2(i),pp2(i)]=ttest2(v1(:,i)-s1(:,i),v2(:,i)-s2(:,i));
    %[pp2(i),hh2(i)]=ranksum(v1(:,i)-s1(:,i),v2(:,i)-s2(:,i));
end
hh2=BenjaminiHochberg(pp2,.1);
hh=BenjaminiHochberg(pp,.1);
fh=figure('Units','Normalized','OuterPosition',[.05 .05 .25 .25]);
subplot(2,1,1)
hold on
p1=plot(mean(v1,1),'LineWidth',2,'DisplayName','VAF');
patch([1:30 30:-1:1],[mean(v1)+std(v1)/sqrt(16) fliplr(mean(v1)-std(v1)/sqrt(16))],p1.Color,'FaceAlpha',.3,'EdgeColor','none')
p2=plot(mean(v2,1),'LineWidth',2,'DisplayName','Patients');
patch([1:30 30:-1:1],[mean(v2)+std(v2)/sqrt(16) fliplr(mean(v2)-std(v2)/sqrt(16))],p2.Color,'FaceAlpha',.3,'EdgeColor','none')
p3=plot(mean(s1,1),'--','LineWidth',2,'Color',p1.Color,'DisplayName','DFT bound');
p4=plot(mean(s2,1),'--','LineWidth',2,'Color',p2.Color,'DisplayName','DFT bound');
plot([1:30], .9*hh, 'k*','MarkerSize',10)
axis([.5 10 .5 .95])
grid on
legend([p1 p3],'Location','SouthEast')
set(gca,'XTickLabel','','FontSize',10,'Position',get(gca,'Position')+[.02 0 0 .03])
ylabel('\bf VAF (\beta)','FontSize',12)
subplot(2,1,2)
hold on
p5=plot(mean(v1-s1,1),'LineWidth',2,'DisplayName','Controls','Color',p1.Color);
p6=plot(mean(v2-s2,1),'LineWidth',2,'DisplayName','Patients','Color',p2.Color);
patch([1:30 30:-1:1],[mean(v1-s1)+std(v1-s1)/sqrt(16) fliplr(mean(v1-s1)-std(v1-s1)/sqrt(16))],p1.Color,'FaceAlpha',.3,'EdgeColor','none')
patch([1:30 30:-1:1],[mean(v2-s2)+std(v2-s2)/sqrt(16) fliplr(mean(v2-s2)-std(v2-s2)/sqrt(16))],p2.Color,'FaceAlpha',.3,'EdgeColor','none')
axis([.5 10 0 .1])
grid on
pp2(pp2>.05)=nan;
plot([1:30], .1*(hh2-.1), 'k*')
aux=mat2cell(num2str([2:2:14]'),ones(7,1),2);
aux{1}='k=2';
set(gca,'FontSize',10,'XTick',[2:2:14],'XTickLabel',aux,'YTick',[0 .05 .1],'YTickLabel',{'0','','0.1'},'Position',get(gca,'Position')+[.02 .03 0 .03])
ylabel('\bf Excess','FontSize',12)
legend([p5 p6],'Location','SouthEast')
%%
saveFig(fh,'./','EMBC2018',0)