%Test merge bounds
%%
bound1=.8-.6./[1:.5:6];
bound1=[.2 .28 .36 .44 .52 .6 .68 .76 .84 .92 .95 .975 .9875 1];
bound2=.6-.3./([1:20].^2);
bound2=[.25 .4 .46 .5:.05:.7 .71 .72 .73 .74 .75 .76 .77 .78 .79 .8];
mixedBound=mergeBounds(bound1,bound2);
%%
figure
subplot(2,1,1)
plot(bound1)
hold on
plot(bound2)
plot(max(bound1,bound2(1:length(bound1))),'LineWidth',2)
plot(mixedBound,'LineWidth',1)
axis([1 length(bound1) 0 1])
grid on

subplot(2,1,2)
plot(diff([0 bound1]))
hold on
plot(diff([0 bound2]))
plot(diff([0 max(bound1,bound2(1:length(bound1)))]),'LineWidth',2)
plot(diff([0 mixedBound]),'LineWidth',1)
axis([1 length(bound1) 0 .4])
grid on
%%