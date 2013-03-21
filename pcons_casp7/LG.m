function P=LG(N,score)
  score=score*20;
  mean=N.^(0.3264)./(0.0437*N.^0.0003+0.0790);
  SD=N./(0.0417*N+3.3700);
  Z=(score/10-mean)./SD;
  expZ=exp(-1*Z);
  if (Z>20)
    P=expZ;
  else 
    P=1-exp(-expZ);
  end
  P=log10(P);