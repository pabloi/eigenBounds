%% See PR to VAF bounds
M=30;

%% Permitted areas:
figure;
subplot(2,1,1)
hold on
cc=get(gca,'ColorOrder');
pr=[2,8,15,25,M];
for i=1:numel(pr)
    [upper,lower]=PR2VAFbound(pr(i),M);
    C=cc(mod(i-1,size(cc,1))+1,:);

    patch([0:M,M:-1:0],[0 upper fliplr(lower) 0],C,'FaceAlpha',.05,'EdgeColor',C,'DisplayName',['PR=' num2str(pr(i))])
end
legend('Location','Southeast')

%% Lower bounds only
subplot(2,1,2)
hold on
cc=get(gca,'ColorOrder');
pr=[1,2,3,5:2:11,15:5:M];
for i=1:numel(pr)
    [upper,lower]=PR2VAFbound(pr(i),M);
    C=cc(mod(i-1,size(cc,1))+1,:);
    plot([0:M],[0 lower],'Color',C,'LineWidth',2,'DisplayName',['PR=' num2str(pr(i))])
end
legend('Location','Southeast')