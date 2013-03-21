N=100:50:1000;
s=100:50:300;

clf
hold on
color=['b','g','r','c','m','k','y'];
for i=1:length(s)
  plot(N,LG(N,s(i)),color(i))
  
end