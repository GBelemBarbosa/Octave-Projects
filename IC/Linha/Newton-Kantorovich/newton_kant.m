A_1 =[13.43538491431723,  -18.14977774439898; 10.01439047867205,  -13.50150142642550]
betam=0.6;
betap=1;
alpham=-0.03305825605413598;

alphap=0.75;
gammam=alpham/betam;
gammap=alphap/betap;
ts = [0.06005004170141785, 0.09007506255212677, 0.225187656380317]
tms = [7.536280233527940, 6.755629691409507, 5.899916597164303]
Ys = [0.791849893551496, 0.820743212919858, 0.928334461670545]
for i=1:3
  Y=Ys(i)
  Ym=Y;
  t=ts(i)
  tm=tms(i)
  A=[0,0,0;0,0,0;0,0,0];
  f=[
  1.533295712619065e-11;
  6.003975094870384e-12;
  3.472999665632415e-12
  ]
  #---------------------------
  A(1,1)=-exp(-alphap*t)*(alphap*Y*sin(t) + (alphap - Y)*cos(t) + sin(t));

  A(1,2)=0;

  A(1,3)=exp(-alphap*t)*sin(t);
  #---------------------------
  A(2,1)=0;

  A(2,2)=exp(alpham*tm)*((A_1(1,1)-A_1(2,2))*cos(betam*tm)/2 + A_1(1,2)*Y*cos(betam*tm) - betam*sin(betam*tm)) + alpham*exp(alpham*tm)*((A_1(1,1)-A_1(2,2))*sin(betam*tm)/(2*betam) + A_1(1,2)*Y*sin(betam*tm)/betam + cos(betam*tm));

  A(2,3)=A_1(1,2)*exp(alpham*tm)*sin(betam*tm)/betam;
  #---------------------------
  A(3,1)=-exp(-alphap*t)*((Y - alphap)*sin(t) + (alphap*Y + 1)*cos(t));

  A(3,2)=-exp(alpham*tm)*(-Y*((A_1(1,1)-A_1(2,2))*cos(betam*tm)/2 + betam*sin(betam*tm)) + A_1(2,1)*cos(betam*tm)) - alpham*exp(alpham*tm)*(Y*(cos(betam*tm) - (A_1(1,1)-A_1(2,2))*sin(betam*tm)/(2*betam)) + A_1(2,1)*sin(betam*tm)/betam);

  A(3,3)=exp(-alphap*t)*cos(t) - exp(alpham*tm)*(cos(betam*tm) - (A_1(1,1)-A_1(2,2))*sin(betam*tm)/(2*betam));
  #---------------------------
  A
  
  B=A\eye(3)
  beta=norm(B, inf)/(1-norm(eye(3)-B*A, inf))
  alpha=beta*f(i)
  #---------------------------
  t-=0.01;
  t=round(t*100)/100;
  t
  tm-=0.01;
  tm=round(tm*100)/100;
  tm
  Y+=0.01;
  Y=round(Y*100)/100;
  Y
  #---------------------------
  exp(-alphap*t)*(abs(alphap^2 - 2*alphap*Y - 1) + abs(Y*alphap^2 + 2*alphap - Y))

  exp(-alphap*t)*(1 + abs(alphap))
  #---------------------------
  abs(2*alpham*exp(alpham*tm)*(abs(A_1(1,2)*Ym) + abs((A_1(1,1)-A_1(2,2))/2) + abs(betam))) + alpham^2*exp(alpham*tm)*(abs(betam*tm) + abs(A_1(1,2)*Ym/betam) + abs((A_1(1,1)-A_1(2,2))/(2*betam))) + exp(alpham*tm)*(betam^2 + abs(betam*A_1(1,2)*Ym) + abs(betam*(A_1(1,1)-A_1(2,2))*sin(betam*tm)/2))

  abs(A_1(1,2)*exp(alpham*tm)*(abs(alpham) + abs(betam))/betam)
  #---------------------------
  exp(-alphap*t)*(abs((-alphap^2 + 2*alphap*Y + 1)) + abs((Y*alphap^2 + 2*alphap - Y)))

  abs(-exp(-alphap*t)*(abs(alphap) + 1))

  #---------------------------
  (abs(2*alpham*exp(alpham*tm))*(abs(A_1(2,1)) + abs((A_1(1,1)-A_1(2,2))*Ym/2) + abs(betam*Ym)) + alpham^2*exp(alpham*tm)*(abs(Ym) + abs(A_1(2,1)/betam) + abs((A_1(1,1)-A_1(2,2))*Ym/(2*betam))) + exp(alpham*tm)*(abs(betam^2*Ym) + abs(betam*A_1(2,1)) + abs(betam*(A_1(1,1)-A_1(2,2))*Ym/2)))

  exp(alpham*tm)*(abs((A_1(2,2)-A_1(1,1))/2) + abs(betam)) + abs(alpham*exp(alpham*tm))*(1 + abs((A_1(2,2)-A_1(1,1))/(2*betam)))
endfor  