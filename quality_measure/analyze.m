A=load('benchmark2/L_4-minsim_3-factor_5.data');
Ssum=A(:,1);
LG=A(:,2);
mx=A(:,3);
time=A(:,4);



A4=load('benchmark2/L_4-minsim_49-factor_2.23.data');
A6=load('benchmark2/L_6-minsim_49-factor_2.23.data');
A8=load('benchmark2/L_6-minsim_49-factor_2.23.data');



Ssum4=A4(:,1);
LG4=A4(:,2);
mx4=A4(:,3);
time4=A4(:,4);

Ssum8=A8(:,1);
LG8=A8(:,2);
mx8=A8(:,3);
time8=A8(:,4);

Ssum6=A6(:,1);
LG6=A6(:,2);
mx6=A6(:,3);
time6=A6(:,4);




LG_all=[LG LG4 LG6 LG8];
Ssum_all=[Ssum Ssum4 Ssum6 Ssum8];
mx_all=[mx mx4 mx6 mx8];
time_all=[time time4 time6 time8];





plot(LG,[LG4,LG6,LG8],'.');

plot(Ssum,[Ssum4,Ssum6+20,Ssum8+40],'.');












