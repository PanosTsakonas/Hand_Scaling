%This code analyses the data extracted from the experimental procedure.
clc;
close all;

D=["Thumb","Index","Middle","Ring","Little"];

%If i am working from my laptop
G=input("are you working from your laptop or uni pc?: ",'s');

if G=='l'

Thumb=importfile("C:\Users\panos\OneDrive - University of Warwick\PhD\Hand Trials\Data for Scaling functions\Hand Scaling.xlsx",D(1));
Index=importfile("C:\Users\panos\OneDrive - University of Warwick\PhD\Hand Trials\Data for Scaling functions\Hand Scaling.xlsx",D(2));
Middle=importfile("C:\Users\panos\OneDrive - University of Warwick\PhD\Hand Trials\Data for Scaling functions\Hand Scaling.xlsx",D(3));
Ring=importfile("C:\Users\panos\OneDrive - University of Warwick\PhD\Hand Trials\Data for Scaling functions\Hand Scaling.xlsx",D(4));
Little=importfile("C:\Users\panos\OneDrive - University of Warwick\PhD\Hand Trials\Data for Scaling functions\Hand Scaling.xlsx",D(5));
M1=importfile("C:\Users\panos\OneDrive - University of Warwick\PhD\Hand Trials\Data for Scaling functions\Hand Scaling.xlsx","Hand Mass");

%If i am working from the University PC
else
Thumb=importfile("C:\Users\u1857308\OneDrive - University of Warwick\PhD\Hand Trials\Data for Scaling functions\Hand Scaling.xlsx",D(1));
Index=importfile("C:\Users\u1857308\OneDrive - University of Warwick\PhD\Hand Trials\Data for Scaling functions\Hand Scaling.xlsx",D(2));
Middle=importfile("C:\Users\u1857308\OneDrive - University of Warwick\PhD\Hand Trials\Data for Scaling functions\Hand Scaling.xlsx",D(3));
Ring=importfile("C:\Users\u1857308\OneDrive - University of Warwick\PhD\Hand Trials\Data for Scaling functions\Hand Scaling.xlsx",D(4));
Little=importfile("C:\Users\u1857308\OneDrive - University of Warwick\PhD\Hand Trials\Data for Scaling functions\Hand Scaling.xlsx",D(5));
M1=importfile("C:\Users\u1857308\OneDrive - University of Warwick\PhD\Hand Trials\Data for Scaling functions\Hand Scaling.xlsx","Hand Mass");
end
n=size(Thumb.data);

%Get the Spring constants and std for each digit
p=1:2:8;
Kt=vpa(Thumb.data(p,1:2:n(2)),5);
SDKt=vpa(Thumb.data(p,2:2:n(2)),5);
Ki=vpa(Index.data(p,1:2:n(2)),5);
SDKi=vpa(Index.data(p,2:2:n(2)),5);
Km=vpa(Middle.data(p,1:2:n(2)),5);
SDKm=vpa(Middle.data(p,2:2:n(2)),5);
Kr=vpa(Ring.data(p,1:2:n(2)),5);
SDKr=vpa(Ring.data(p,2:2:n(2)),5);
Kl=vpa(Little.data(p,1:2:n(2)),5);
SDKl=vpa(Little.data(p,2:2:n(2)),5);

%Get the damper constants and std for each digit
p1=2:2:8;
Bt=vpa(Thumb.data(p1,1:2:n(2)),5);
SDBt=vpa(Thumb.data(p1,2:2:n(2)),5);
Bi=vpa(Index.data(p1,1:2:n(2)),5);
SDBi=vpa(Index.data(p1,2:2:n(2)),5);
Bm=vpa(Middle.data(p1,1:2:n(2)),5);
SDBm=vpa(Middle.data(p1,2:2:n(2)),5);
Br=vpa(Ring.data(p1,1:2:n(2)),5);
SDBr=vpa(Ring.data(p1,2:2:n(2)),5);
Bl=vpa(Little.data(p1,1:2:n(2)),5);
SDBl=vpa(Little.data(p1,2:2:n(2)),5);

%Get the lengths and radii from the data 
%The 1 collumn corresponds to MCP 2 Pip and 3 DIP
Lt1=Thumb.data(11:n(1),1:2:6);
Lt=[Lt1(:,3) Lt1(:,2) Lt1(:,1)];
Rt1=Thumb.data(11:n(1),2:2:6);
Rt=[Rt1(:,3) Rt1(:,2) Rt1(:,1)];
Mt=Lt.*Rt.^2.*(pi*1.1*10^-3);
M=Mt.*10^-3;
Icmc=M(:,1).*(Rt(:,1).^2./4+Lt(:,1).^2./3).*10^-6+(M(:,2)+M(:,3)).*Lt(:,1).^2.*10^-6;
Imcp=M(:,2).*(Rt(:,2).^2./4+Lt(:,2).^2./3).*10^-6+M(:,3).*Lt(:,2).^2.*10^-6;
Iip=M(:,3).*(Rt(:,3).^2./4+Lt(:,3).^2./3).*10^-6;
clear M;
Li1=Index.data(11:n(1),1:2:6);
Li=[Li1(:,3) Li1(:,2) Li1(:,1)];
Ri1=Index.data(11:n(1),2:2:6);
Ri=[Ri1(:,3) Ri1(:,2) Ri1(:,1)];
Mi=Li.*Ri.^2.*(pi*1.1*10^-3);
M=Mi.*10^-3;
Imcpi=M(:,1).*(Ri(:,1).^2./4+Li(:,1).^2./3).*10^-6+(M(:,2)+M(:,3)).*Li(:,1).^2.*10^-6;
Ipipi=M(:,2).*(Ri(:,2).^2./4+Li(:,2).^2./3).*10^-6+M(:,3).*Li(:,2).^2.*10^-6;
Idipi=M(:,3).*(Ri(:,3).^2./4+Li(:,3).^2./3).*10^-6;
clear M;
Lm1=Middle.data(11:n(1),1:2:6);
Lm=[Lm1(:,3) Lm1(:,2) Lm1(:,1)];
Rm1=Middle.data(11:n(1),2:2:6);
Rm=[Rm1(:,3) Rm1(:,2) Rm1(:,1)];
Mm=Lm.*Rm.^2.*(pi*1.1*10^-3);
M=Mm.*10^-3;
Imcpm=M(:,1).*(Rm(:,1).^2./4+Lm(:,1).^2./3).*10^-6+(M(:,2)+M(:,3)).*Lm(:,1).^2.*10^-6;
Ipipm=M(:,2).*(Rm(:,2).^2./4+Lm(:,2).^2./3).*10^-6+M(:,3).*Lm(:,2).^2.*10^-6;
Idipm=M(:,3).*(Rm(:,3).^2./4+Lm(:,3).^2./3).*10^-6;
clear M;
Lr1=Ring.data(11:n(1),1:2:6);
Lr=[Lr1(:,3) Lr1(:,2) Lr1(:,1)];
Rr1=Ring.data(11:n(1),2:2:6);
Rr=[Rr1(:,3) Rr1(:,2) Rr1(:,1)];
Mr=Lr.*Rr.^2.*(pi*1.1*10^-3);
M=Mr.*10^-3;
Imcpr=M(:,1).*(Rr(:,1).^2./4+Lr(:,1).^2./3).*10^-6+(M(:,2)+M(:,3)).*Lr(:,1).^2.*10^-6;
Ipipr=M(:,2).*(Rr(:,2).^2./4+Lr(:,2).^2./3).*10^-6+M(:,3).*Lr(:,2).^2.*10^-6;
Idipr=M(:,3).*(Rr(:,3).^2./4+Lr(:,3).^2./3).*10^-6;
clear M;
Ll1=Little.data(11:n(1),1:2:6);
Ll=[Ll1(:,3) Ll1(:,2) Ll1(:,1)];
Rl1=Little.data(11:n(1),2:2:6);
Rl=[Rl1(:,3) Rl1(:,2) Rl1(:,1)];
Ml=Ll.*Rl.^2.*(pi*1.1*10^-3);
M=Ml.*10^-3;
Imcpl=M(:,1).*(Rl(:,1).^2./4+Ll(:,1).^2./3).*10^-6+(M(:,2)+M(:,3)).*Ll(:,1).^2.*10^-6;
Ipipl=M(:,2).*(Rl(:,2).^2./4+Ll(:,2).^2./3).*10^-6+M(:,3).*Ll(:,2).^2.*10^-6;
Idipl=M(:,3).*(Rl(:,3).^2./4+Ll(:,3).^2./3).*10^-6;
clear M;
clear Lt1 Rt1 Li1 Ri1 Lm1 Rm1 Lr1 Rr1 Ll1 Rl1
I=[Icmc Imcp Iip Imcpi Ipipi Idipi Imcpm Ipipm Idipm Imcpr Ipipr Idipr Imcpl Ipipl Idipl];

%Palm Length, Hand breadth and Hand width
PL=Thumb.data(11:n(1),7);
HB=Thumb.data(11:n(1),8);
HW=Thumb.data(11:n(1),9);
HL=PL+Lm(:,1)+Lm(:,2)+Lm(:,3);
Palm=PL+HB+HW;
PLHB=PL.*HB;
HLHB=HL.*HB;

%Get the theoretical and experimental hand masses
M_th=M1.data(:,1);
M_exp=M1.data(:,2);

%Anthropometric tables
g=Lt.*Rt;
g1=Li.*Ri;
g2=Lm.*Rm;
g3=Lr.*Rr;
g4=Ll.*Rl;
L=[Lt Li Lm Lr Ll ];
La=[g(:,1)+g(:,2)+g(:,3) g1(:,1)+g1(:,2)+g1(:,3) g2(:,1)+g2(:,2)+g2(:,3) g3(:,1)+g3(:,2)+g3(:,3) g4(:,1)+g4(:,2)+g4(:,3)];
Ra=[Rt Ri Rm Rr Rl];


%The next function calculates the Bland Altman plot and it checks if the mean difference between the 2 methods
%follows a gaussian distribution using the Kolmogorov-Smirnov test. If the data seems to violate the assumption of 
%distribution type, a warning message is generated.

BlandAltman(M_exp,M_th,{'Experimental','Theoretical', 'Grams'},'corrInfo',{'eq';'r2';'r';'n';'RMSE'},'baInfo',{'LOA';'ks'},'axesLimits','tight','baStatsMode','NORMAL', 'showFitCI','off');

%Correlation coefficients
[rPear,p_v,RL,RU]=corrcoef(M_th(1:23),M_exp(1:23),'alpha',0.05);

fitfun=fittype(@(a,b,x) a+b.*x);
GoF_RMSE=zeros(15,1);
x0=[0.001 0.2031];
SL=zeros(23,15);
for j=1:15
    
[Fit,GoF]=fit(Palm,L(:,j),fitfun,'MaxFunEvals',10^8,'MaxIter',10^8,'ToLFun',10^-8,'Robust','Bisquare','StartPoint',x0);
[F,G]=fit(PL,L(:,j),fitfun,'MaxFunEvals',10^8,'MaxIter',10^8,'ToLFun',10^-8,'Robust','Bisquare','StartPoint',x0);
[F1,G1]=fit(PL.*HB,L(:,j),fitfun,'MaxFunEvals',10^8,'MaxIter',10^8,'ToLFun',10^-8,'Robust','Bisquare','StartPoint',x0);
[F2,G2]=fit(HL,L(:,j),fitfun,'MaxFunEvals',10^8,'MaxIter',10^8,'ToLFun',10^-8,'Robust','Bisquare','StartPoint',x0);
[F3,G3]=fit(HB,L(:,j),fitfun,'MaxFunEvals',10^8,'MaxIter',10^8,'ToLFun',10^-8,'Robust','Bisquare','StartPoint',x0);
[F4,G4]=fit(HW,L(:,j),fitfun,'MaxFunEvals',10^8,'MaxIter',10^8,'ToLFun',10^-8,'Robust','Bisquare','StartPoint',x0);
LM=fitlm([HL HB HW],L(:,j));
LM1=fitlm([PL HB HW],L(:,j));

[Fit1,GoF1]=fit(Palm,Ra(:,j),fitfun,'MaxFunEvals',10^8,'MaxIter',10^8,'ToLFun',10^-8,'Robust','Bisquare','StartPoint',x0);
[Ft1,GOFt1]=fit(HB,Ra(:,j),fitfun,'MaxFunEvals',10^8,'MaxIter',10^8,'ToLFun',10^-8,'Robust','Bisquare','StartPoint',x0);
[Ft2,GOFt2]=fit(PL.*HB,Ra(:,j),fitfun,'MaxFunEvals',10^8,'MaxIter',10^8,'ToLFun',10^-8,'Robust','Bisquare','StartPoint',x0);
[Ft3,GOFt3]=fit(HL,Ra(:,j),fitfun,'MaxFunEvals',10^8,'MaxIter',10^8,'ToLFun',10^-8,'Robust','Bisquare','StartPoint',x0);
[Ft4,GOFt4]=fit(HW,Ra(:,j),fitfun,'MaxFunEvals',10^8,'MaxIter',10^8,'ToLFun',10^-8,'Robust','Bisquare','StartPoint',x0);
[Ft5,GOFt5]=fit(PL,Ra(:,j),fitfun,'MaxFunEvals',10^8,'MaxIter',10^8,'ToLFun',10^-8,'Robust','Bisquare','StartPoint',x0);
RM=fitlm([HL HB HW],Ra(:,j));
RM1=fitlm([PL HB HW],Ra(:,j));

GOF_L=max([GoF.rsquare G.rsquare G1.rsquare G2.rsquare G3.rsquare G4.rsquare LM.Rsquared.Ordinary LM1.Rsquared.Ordinary]);
GG=[GoF.rsquare G.rsquare G1.rsquare G3.rsquare G4.rsquare LM1.Rsquared.Ordinary];
GOF_R=max([GoF1.rsquare GOFt1.rsquare GOFt2.rsquare GOFt3.rsquare GOFt4.rsquare GOFt5.rsquare RM.Rsquared.Ordinary RM1.Rsquared.Ordinary]);

if j<=3
    f="thumb";
elseif j<=6
    f="index";
elseif j<=9
    f="middle";
elseif j<=12
    f="ring";
else
    f="little";
end

if mod(j,3)==1 && j<=3
    g="metacarpal";
else
    g="proximal";
end
if mod(j,3)==2
    g="middle";
elseif mod(j,3)==0
    g="distal";
end
h=figure('Position',get(0,'Screensize'));


if j>6 && j<=9
    
subplot(2,1,1)
if max(GG)==(GoF.rsquare)
    coef=coeffvalues(Fit);
    SL(:,j)=Fit(Palm);
plot(Palm,L(:,j),'o',Palm,Fit(Palm));
c="y= "+num2str(coef(2))+"*x+"+num2str(coef(1));
title("Linear regression for the length of "+f+" "+g+" segment with R^2= "+GoF.rsquare);
xlabel("Sum of PL, HB and HW dimensions (mm)");
ylabel("Segment length (mm)");
legend(c);


elseif max(GG)==(G.rsquare)
    coef=coeffvalues(F);
    SL(:,j)=F(PL);
plot(PL,L(:,j),'o',PL,F(PL));
c="y= "+num2str(coef(2))+"*x+"+num2str(coef(1))+" with RMSE: "+G.rmse+" mm";
title("Linear regression for the length of "+f+" "+g+" segment with R^2= "+G.rsquare);
xlabel("Palm Length (mm)");
ylabel("Segment length (mm)");
legend(c);
elseif max(GG)==(G1.rsquare)
   coef=coeffvalues(F1);
   SL(:,j)=F1(PL.*HB);
plot(PL.*HB,L(:,j),'o',PL.*HB,F1(PL.*HB));
c="y= "+num2str(coef(2))+"*x+"+num2str(coef(1))+" with RMSE: "+G1.rmse+" mm";
title("Linear regression for the length of "+f+" "+g+" segment with R^2= "+G1.rsquare);
xlabel("Product of Palm Length and Hand breadth (mm^2)");
ylabel("Segment length (mm)");
legend(c); 
elseif max(GG)==(G3.rsquare)
   coef=coeffvalues(F3);
   SL(:,j)=F3(HB);
plot(HB,L(:,j),'o',HB,F3(HB));
c="y= "+num2str(coef(2))+"*x+"+num2str(coef(1))+" with RMSE: "+G3.rmse+" mm";
title("Linear regression for the length of "+f+" "+g+" segment with R^2= "+G3.rsquare);
xlabel("Hand breadth (mm)");
ylabel("Segment length (mm)");
legend(c);
elseif max(GG)==(G4.rsquare)
    coef=coeffvalues(F4);
    SL(:,j)=F4(HW);
plot(HW,L(:,j),'o',HW,F4(HW));
c="y= "+num2str(coef(2))+"*x+"+num2str(coef(1))+" with RMSE: "+G4.rmse+" mm";
title("Linear regression for the length of "+f+" "+g+" segment with R^2= "+G4.rsquare);
xlabel("Hand width (mm)");
ylabel("Segment length (mm)");
legend(c);
else
    plot3(PL,HB,L(:,j),'o');
    hold on;
    coef=table2array(LM1.Coefficients(:,1));
    SL(:,j)=coef(1)+coef(2).*PL+coef(3).*HB+coef(4).*HW;
    x1f=linspace(min(PL),max(PL),20);
    x2f=linspace(min(HB),max(HB),20);
    [X1F,X2F]=meshgrid(x1f,x2f);
    Yf=coef(1)+coef(2)*X1F+coef(3)*X2F;
    mesh(X1F,X2F,Yf)
    hold off;
    xlabel("Palm Length (mm)");
    ylabel("Hand Breadth (mm)");
    zlabel("Segment Length (mm)");
    c="y= "+num2str(coef(1))+"+PL*"+num2str(coef(2))+"+HB*"+num2str(coef(3))+"+HW*"+num2str(coef(4))+" with RMSE: "+LM1.RMSE+" mm";
    title("Multivariate Linear regression for the length of "+f+" "+g+" segment with R^2= "+LM1.Rsquared.Ordinary);
    legend(c);
    
end
else
    subplot(2,1,1)
if GOF_L==(GoF.rsquare)
    coef=coeffvalues(Fit);
    SL(:,j)=Fit(Palm);
plot(Palm,L(:,j),'o',Palm,Fit(Palm));
c="y= "+num2str(coef(2))+"*x+"+num2str(coef(1))+" with RMSE: "+GoF.rmse+" mm";
title("Linear regression for the length of "+f+" "+g+" segment with R^2= "+GoF.rsquare);
xlabel("Sum of PL, HB and HW dimensions (mm)");
ylabel("Segment length (mm)");
legend(c);

elseif GOF_L==(G.rsquare)
    coef=coeffvalues(F);
    SL(:,j)=F(PL);
plot(PL,L(:,j),'o',PL,F(PL));
c="y= "+num2str(coef(2))+"*x+"+num2str(coef(1))+" with RMSE: "+G.rmse+" mm";
title("Linear regression for the length of "+f+" "+g+" segment with R^2= "+G.rsquare);
xlabel("Palm Length (mm)");
ylabel("Segment length (mm)");
legend(c);
elseif GOF_L==(G1.rsquare)
   coef=coeffvalues(F1);
   SL(:,j)=F1(PL.*HB);
plot(PL.*HB,L(:,j),'o',PL.*HB,F1(PL.*HB));
c="y= "+num2str(coef(2))+"*x+"+num2str(coef(1))+" with RMSE: "+G1.rmse+" mm";
title("Linear regression for the length of "+f+" "+g+" segment with R^2= "+G1.rsquare);
xlabel("Product of Palm Length and Hand breadth (mm^2)");
ylabel("Segment length (mm)");
legend(c); 
elseif GOF_L==(G2.rsquare)
coef=coeffvalues(F2);
SL(:,j)=F2(HL);
plot(HL,L(:,j),'o',HL,F2(HL));
c="y= "+num2str(coef(2))+"*x+"+num2str(coef(1))+" with RMSE: "+G2.rmse+" mm";
title("Linear regression for the length of "+f+" "+g+" segment with R^2= "+G2.rsquare);
xlabel("Hand length (mm)");
ylabel("Segment length (mm)");
legend(c);     
elseif GOF_L==(G3.rsquare)
coef=coeffvalues(F3);
SL(:,j)=F3(HB);
plot(HB,L(:,j),'o',HB,F3(HB));
c="y= "+num2str(coef(2))+"*x+"+num2str(coef(1))+" with RMSE: "+G3.rmse+" mm";
title("Linear regression for the length of "+f+" "+g+" segment with R^2= "+G3.rsquare);
xlabel("Hand breadth (mm)");
ylabel("Segment length (mm)");
legend(c); 
elseif GOF_L==(G4.rsquare)
    coef=coeffvalues(F4);
    SL(:,j)=F4(HW);
plot(HW,L(:,j),'o',HW,F4(HW));
c="y= "+num2str(coef(2))+"*x+"+num2str(coef(1))+" with RMSE: "+G4.rmse+" mm";
title("Linear regression for the length of "+f+" "+g+" segment with R^2= "+G4.rsquare);
xlabel("Hand width (mm)");
ylabel("Segment length (mm)");
legend(c);
elseif GOF_L==LM.Rsquared.Ordinary
    plot3(HL,HB,L(:,j),'o');
    hold on;
    coef=table2array(LM.Coefficients(:,1));
    SL(:,j)=coef(1)+coef(2).*HL+coef(3).*HB+coef(4).*HW;
    x1f=linspace(min(HL),max(HL),20);
    x2f=linspace(min(HB),max(HB),20);
    [X1F,X2F]=meshgrid(x1f,x2f);
    Yf=coef(1)+coef(2)*X1F+coef(3)*X2F;
    mesh(X1F,X2F,Yf)
    xlabel("Hand Length (mm)");
    ylabel("Hand Breadth (mm)");
    zlabel("Segment Length (mm)");
    c="y= "+num2str(coef(1))+"+HL*"+num2str(coef(2))+"+HB*"+num2str(coef(3))+"+HW*"+num2str(coef(4))+" with RMSE: "+LM.RMSE+" mm";
    title("Multivariate Linear regression for the length of "+f+" "+g+" segment with R^2= "+LM.Rsquared.Ordinary);
    legend(c);
else
    plot3(PL,HB,L(:,j),'o');
    hold on;
    coef=table2array(LM1.Coefficients(:,1));
    SL(:,j)=coef(1)+coef(2).*PL+coef(3).*HB+coef(4).*HW;
    x1f=linspace(min(PL),max(PL),20);
    x2f=linspace(min(HB),max(HB),20);
    [X1F,X2F]=meshgrid(x1f,x2f);
    Yf=coef(1)+coef(2)*X1F+coef(3)*X2F;
    mesh(X1F,X2F,Yf)
    hold off;
    xlabel("Palm Length (mm)");
    ylabel("Hand Breadth (mm)");
    zlabel("Segment Length (mm)");
    c="y= "+num2str(coef(1))+"+PL*"+num2str(coef(2))+"+HB*"+num2str(coef(3))+"+HW*"+num2str(coef(4))+" with RMSE: "+LM1.RMSE+" mm";
    title("Multivariate Linear regression for the length of "+f+" "+g+" segment with R^2= "+LM1.Rsquared.Ordinary);
    legend(c);
end
end

subplot(2,1,2)
if GOF_R==(GoF1.rsquare)
    coef1=coeffvalues(Fit1);
plot(Palm,Ra(:,j),'x',Palm,Fit1(Palm));
ylabel("Segment radius (mm)");
xlabel("Sum of PL, HB and HW dimensions (mm)");
c1="y= "+num2str(coef1(2))+"*x+"+num2str(coef1(1))+" with RMSE: "+GoF1.rmse+" mm";
title("Linear regression for the radius of "+f+" "+g+" segment with R^2= "+GoF1.rsquare);
legend(c1);
elseif GOF_R==(GOFt1.rsquare)
    coef1=coeffvalues(Ft1);
  plot(HB,Ra(:,j),'x',HB,Ft1(HB));
ylabel("Segment radius (mm)");
xlabel("Hand breadth (mm)");
c1="y= "+num2str(coef1(2))+"*x+"+num2str(coef1(1))+" with RMSE: "+GOFt1.rmse+" mm";
title("Linear regression for the radius of "+f+" "+g+" segment with R^2= "+GOFt1.rsquare);
legend(c1);  
elseif GOF_R==(GOFt2.rsquare)
    coef1=coeffvalues(Ft2);
  plot(PL.*HB,Ra(:,j),'x',PL.*HB,Ft2(PL.*HB));
ylabel("Segment radius (mm)");
xlabel("Product of palm length and hand breadth (mm^2)");
c1="y= "+num2str(coef1(2))+"*x+"+num2str(coef1(1))+" with RMSE: "+GOFt2.rmse+" mm";
title("Linear regression for the radius of "+f+" "+g+" segment with R^2= "+GOFt2.rsquare);
legend(c1);  
elseif GOF_R==(GOFt3.rsquare)
    coef1=coeffvalues(Ft3);
  plot(HL,Ra(:,j),'x',HL,Ft3(HL));
ylabel("Segment radius (mm)");
xlabel("Hand length (mm)");
c1="y= "+num2str(coef1(2))+"*x+"+num2str(coef1(1))+" with RMSE: "+GOFt3.rmse+" mm";
title("Linear regression for the radius of "+f+" "+g+" segment with R^2= "+GOFt3.rsquare);
legend(c1); 
elseif GOF_R==(GOFt4.rsquare)
    coef1=coeffvalues(Ft4);
  plot(HW,Ra(:,j),'x',HW,Ft4(HW));
ylabel("Segment radius (mm)");
xlabel("Hand width (mm)");
c1="y= "+num2str(coef1(2))+"*x+"+num2str(coef1(1))+" with RMSE: "+GOFt4.rmse+" mm";
title("Linear regression for the radius of "+f+" "+g+" segment with R^2= "+GOFt4.rsquare);
legend(c1); 
elseif GOF_R==(GOFt5.rsquare)
    coef1=coeffvalues(Ft5);
  plot(PL,Ra(:,j),'x',PL,Ft5(PL));
ylabel("Segment radius (mm)");
xlabel("Palm length (mm)");
c1="y= "+num2str(coef1(2))+"*x+"+num2str(coef1(1))+" with RMSE: "+GOFt5.rmse+" mm";
title("Linear regression for the radius of "+f+" "+g+" segment with R^2= "+GOFt5.rsquare);
legend(c1); 
elseif GOF_R==RM.Rsquared.Ordinary
    plot3(HL,HB,Ra(:,j),'o');
    hold on;
    coef=table2array(RM.Coefficients(:,1));
    x1f=linspace(min(HL),max(HL),20);
    x2f=linspace(min(HB),max(HB),20);
    [X1F,X2F]=meshgrid(x1f,x2f);
    Yf=coef(1)+coef(2)*X1F+coef(3)*X2F;
    mesh(X1F,X2F,Yf)
    hold off;
    xlabel("Hand Length (mm)");
    ylabel("Hand Breadth (mm)");
    zlabel("Segment Radius (mm)");
    c="y= "+num2str(coef(1))+"+HL*"+num2str(coef(2))+"+HB*"+num2str(coef(3))+"+HW*"+num2str(coef(4))+" with RMSE: "+RM.RMSE+" mm";
    title("Multivariate Linear regression for the radius of "+f+" "+g+" segment with R^2= "+RM.Rsquared.Ordinary);
    legend(c);
else
    plot3(PL,HB,Ra(:,j),'o');
    hold on;
    coef=table2array(RM1.Coefficients(:,1));
    x1f=linspace(min(PL),max(PL),20);
    x2f=linspace(min(HB),max(HB),20);
    [X1F,X2F]=meshgrid(x1f,x2f);
    Yf=coef(1)+coef(2)*X1F+coef(3)*X2F;
    mesh(X1F,X2F,Yf)
    hold off;
    xlabel("Palm Length (mm)");
    ylabel("Hand Breadth (mm)");
    zlabel("Segment Radius (mm)");
    c="y= "+num2str(coef(1))+"+PL*"+num2str(coef(2))+"+HB*"+num2str(coef(3))+"+HW*"+num2str(coef(4))+" with RMSE: "+RM1.RMSE+" mm";
    title("Multivariate Linear regression for the radius of "+f+" "+g+" segment with R^2= "+RM1.Rsquared.Ordinary);
    legend(c);
end
saveas(h,f+" "+g,'jpg');
clear LM LM1 RM RM1 Fit F1 F2 F3 F4 Fit1 Ft1 Ft2 Ft3 Ft4 Ft5;
end
%}

%This part calculates the difference between estimating the segment length
%using the scaling functions determined here and the ones from Buchholz

syms x
LMtB=feval(matlabFunction(0.251*x),HL);
LPtB=feval(matlabFunction(0.196*x),HL);
LDtB=feval(matlabFunction(0.158*x),HL);
LPiB=feval(matlabFunction(0.245*x),HL);
LMiB=feval(matlabFunction(0.143*x),HL);
LDiB=feval(matlabFunction(0.097*x),HL);
LPmB=feval(matlabFunction(0.266*x),HL);
LMmB=feval(matlabFunction(0.170*x),HL);
LDmB=feval(matlabFunction(0.108*x),HL);
LPrB=feval(matlabFunction(0.244*x),HL);
LMrB=feval(matlabFunction(0.165*x),HL);
LDrB=feval(matlabFunction(0.107*x),HL);
LPlB=feval(matlabFunction(0.204*x),HL);
LMlB=feval(matlabFunction(0.117*x),HL);
LDlB=feval(matlabFunction(0.093*x),HL);

BB=[LMtB LPtB LDtB LPiB LMiB LDiB LPmB LMmB LDmB LPrB LMrB LDrB LPlB LMlB LDlB]; 
writematrix(BB,'Buchholz.csv');
writematrix(SL,'This_study.csv');
DB=zeros(23,15);
DSL=DB;
for j=1:15
    
    if j<=3
    f="thumb";
elseif j<=6
    f="index";
elseif j<=9
    f="middle";
elseif j<=12
    f="ring";
else
    f="little";
    end

if mod(j,3)==1 && j<=3
    g="metacarpal";
else
    g="proximal";
end
if mod(j,3)==2
    g="middle";
elseif mod(j,3)==0
    g="distal";
end
    h=figure();
    
        mm=stem(L(:,j),BB(:,j),'o','-.');
        mm.LineWidth=2;
        grid
        hold on;
        mm1=stem(L(:,j),SL(:,j),'o',':');
        mm1.LineWidth=2;
        hold on;
       grid
        xlabel("Participant segment length data (mm)");
        ylabel("Predicted segment length (mm)");
        title({'Graph between the predicted segment length values of '+f+' finger',...
                ' '+g+' segment from this study, the Buchholz scaling function and the participant data'});
        legend("Buchholz function","Scaling function from this study",'Location','southeast');
        ax=gca;
        ax.FontSize=8;
    hold off;
    saveas(h,"Predicted vs participant "+f+" "+g,'jpg');
    saveas(h,"Predicted vs participant "+f+" "+g,'fig');
    clear h;
    h=figure();
    DB(:,j)=L(:,j)-BB(:,j);
    DSL(:,j)=L(:,j)-SL(:,j);
    plot(L(:,j),DB(:,j),'o',L(:,j),DSL(:,j),'x');
    xlabel("Participant segment length data (mm)");
    ylabel("Difference between participant data and predicted value (mm)");
    legend("Buchholz equation","Scaling equation from this study",'Location','southeast');
    title({'Difference between participant data and scaling equation',... 
            'for '+f+' finger '+g+' segment'});
     ax=gca;
     ax.FontSize=8;
    saveas(h,"Scaling equation Difference "+f+" "+g,'jpg');
    saveas(h,"Scaling equation Difference "+f+" "+g,'fig');
end


%Non Linear model fitting
model= @(a,x) a(1)+a(2).*x;
R=zeros(length(Li),2);
opts = statset('TolFun',1e-10,'MaxIter',10^10);
%Set the standard deviation and the spring and damper constnats for the
%flexion/extension
SDK=[SDKt(1:3,:);SDKi(1:3,:);SDKm(1:3,:);SDKr(1:3,:);SDKl(1:3,:)];
SDB=[SDBt(1:3,:);SDBi(1:3,:);SDBm(1:3,:);SDBr(1:3,:);SDBl(1:3,:)];
K=[Kt(1:3,:);Ki(1:3,:);Km(1:3,:);Kr(1:3,:);Kl(1:3,:)];
B=[Bt(1:3,:);Bi(1:3,:);Bm(1:3,:);Br(1:3,:);Bl(1:3,:)];
%Set the standard deviation and the spring and damper constnats for the
%ab/adduction
Ka=[Kt(4,:);Ki(4,:);Km(4,:);Kr(4,:);Kl(4,:)];
Ba=[Bt(4,:);Bi(4,:);Bm(4,:);Br(4,:);Bl(4,:)];
SDKa=[SDKt(4,:);SDKi(4,:);SDKm(4,:);SDKr(4,:);SDKl(4,:)];
SDBa=[SDBt(4,:);SDBi(4,:);SDBm(4,:);SDBr(4,:);SDBl(4,:)];

%Calculates mean damping ratio for digits in flexion/extension
z=zeros(2,15);
for i=1:15
    z(1,i)=mean(double(B(i,:))./(2.*sqrt(double(K(i,:)).*transpose(I(:,i)))),'omitnan');
    z(2,i)=std(double(B(i,:))./(2.*sqrt(double(K(i,:)).*transpose(I(:,i)))),'omitnan');
end

%Calculates mean damping ration for abduction movement

z1=zeros(2,5);
Iab=zeros(5,23);
for i=1:5
    for j=1:23
        %From Laptop
        %Ia=importfile("C:\Users\panos\OneDrive - University of Warwick\PhD\Hand Trials\Digit Weights\P_"+j+".xlsx","");
        %From computer at the University
       Ia=importfile("C:\Users\u1857308\OneDrive - University of Warwick\PhD\Hand Trials\Digit Weights\P_"+j+".xlsx","");
        Iab(i,j)=Ia.data(27+i,7);
        clear Ia;
    end
    z1(1,i)=mean(double(Ba(i,:))./(2.*sqrt(double(Ka(i,:)).*Iab(i,:))),'omitnan');
    z1(2,i)=std(double(Ba(i,:))./(2.*sqrt(double(Ka(i,:)).*Iab(i,:))),'omitnan');
end

J=1:3:15;
J1=2:3:15;
J2=3:3:15;

if isempty(j)==1
    j=1;
end

%Flexion-Extension scaling functions
for i=1:15
%Add the weights as the reciprocal of the standard deviation
w=1./SDK(i,:);
wb=1./SDB(i,:);
wb=double(wb./max(wb));
wk=double(w./max(w));

    
ff=fitnlm(L(:,i).*Ra(:,i),double(K(i,:)),model,[0.003 0.01],'Weight',transpose(wk),'Options',opts);
f1=fitnlm(L(:,i),double(K(i,:)),model,[0.003 0.01],'Weight',transpose(wk),'Options',opts);
f2=fitnlm(PL,double(K(i,:)),model,[0.003 0.01],'Weight',transpose(wk),'Options',opts);
f3=fitnlm(PL.*HB,double(K(i,:)),model,[0.003 0.01],'Weight',transpose(wk),'Options',opts);
f4=fitnlm(HLHB,double(K(i,:)),model,[0.003 0.01],'Weight',transpose(wk),'Options',opts);
f5=fitnlm(Ra(:,i),double(K(i,:)),model,[0.003 0.01],'Weight',transpose(wk),'Options',opts);
f6=fitnlm(Palm,double(K(i,:)),model,[0.003 0.01],'Weight',transpose(wk),'Options',opts);

ff1=fitnlm(L(:,i).*Ra(:,i),double(B(i,:)),model,[0.003 0.01],'Weight',transpose(wb),'Options',opts);
ff2=fitnlm(Ra(:,i),double(B(i,:)),model,[0.003 0.01],'Weight',transpose(wb),'Options',opts);
ff3=fitnlm(PL,double(B(i,:)),model,[0.003 0.01],'Weight',transpose(wb),'Options',opts);
ff4=fitnlm(PL.*HB,double(B(i,:)),model,[0.003 0.01],'Weight',transpose(wb),'Options',opts);
ff5=fitnlm(HLHB,double(B(i,:)),model,[0.003 0.01],'Weight',transpose(wb),'Options',opts);
ff6=fitnlm(L(:,i),double(B(i,:)),model,[0.003 0.01],'Weight',transpose(wb),'Options',opts);
ff7=fitnlm(Palm,double(B(i,:)),model,[0.003 0.01],'Weight',transpose(wb),'Options',opts);
  R(i,1)=ff.Rsquared.Ordinary;
  R(i,2)=ff1.Rsquared.Ordinary;

[Rf,kl]=max([ff.Rsquared.Ordinary f1.Rsquared.Ordinary f2.Rsquared.Ordinary f3.Rsquared.Ordinary f4.Rsquared.Ordinary f5.Rsquared.Ordinary f6.Rsquared.Ordinary]);
[Rff,kl1]=max([ff1.Rsquared.Ordinary ff2.Rsquared.Ordinary ff3.Rsquared.Ordinary ff4.Rsquared.Ordinary ff5.Rsquared.Ordinary ff6.Rsquared.Ordinary ff7.Rsquared.Ordinary]);
 Rms=[ff.RMSE f1.RMSE f2.RMSE f3.RMSE f4.RMSE f5.RMSE f6.RMSE];
 Rms1=[ff1.RMSE ff2.RMSE ff3.RMSE ff4.RMSE ff5.RMSE ff6.RMSE ff7.RMSE];
b11=ispart(i,J);
b1=ispart(i,J1);
b2=ispart(i,J2);

h=figure('Position',get(0,'Screensize'));

if i<=3
    p="thumb";
elseif i<=6
    p="index";
elseif i<=9
    p="middle";
    elseif i<=12
        p="ring";
else
    p="little";
end

if b11==1
    ck="Spring regression for MCP joint of "+p+" finger with R^2 "+Rf+" and RMSE "+Rms(kl)+" Nm/rad";
    cb="Damper regression for MCP joint of "+p+" finger with R^2 "+Rff+" and RMSE "+Rms1(kl1)+" Nms/rad";
    ShName="MCP "+p;
    cp="Non linear model fit for MCP joint of "+p;
    lk="Spring MCP data of "+p;
    ld="Damper MCP data of "+p;
elseif b1==1
    ck="Spring regression for PIP joint of "+p+" finger with R^2 "+Rf+" and RMSE "+Rms(kl)+" Nm/rad";
    cb="Damper regression for PIP joint of "+p+" finger with R^2 "+Rff+" and RMSE "+Rms1(kl1)+" Nms/rad";
    ShName="PIP "+p;
    cp="Non linear model fit for PIP joint of "+p;
    lk="Spring PIP data of "+p;
    ld="Damper PIP data of "+p;
elseif b2==1
    ck="Spring regression for DIP joint of "+p+" finger with R^2 "+Rf+" and RMSE "+Rms(kl)+" Nm/rad";
    cb="Damper regression for DIP joint of "+p+" finger with R^2 "+Rff+" and RMSE "+Rms1(kl1)+" Nms/rad";
    ShName="DIP "+p;
    cp="Non linear model fit for DIP joint of "+p;
    lk="Spring DIP data of "+p;
    ld="Damper DIP data of "+p;
end

if Rf== ff.Rsquared.Ordinary 
subplot(2,1,1)
errorbar(L(:,i).*Ra(:,i),K(i,:),SDK(i,:),'o');
hold on;
if Rf<0.5
plot(L(:,i).*Ra(:,i),predict(ff,L(:,i).*Ra(:,i)),L(:,i).*Ra(:,i),mean(double(K(i,:)),'omitnan').*ones(1,length(PL)),'--');
legend(lk,"Fit "+ff.Coefficients{1,1}+" + "+ff.Coefficients{2,1}+" *x", "Mean Spring value "+mean(double(K(i,:)),'omitnan'),'Location','southeast');
else
    plot(L(:,i).*Ra(:,i),predict(ff,L(:,i).*Ra(:,i)));
    legend(lk,"Fit "+ff.Coefficients{1,1}+" + "+ff.Coefficients{2,1}+" *x",'Location','southeast');
end

ylabel("K (Nm/rad)");
xlabel("Product of segment length and radius (mm^2)");
title(ck);
hold off;

writematrix("y=a1+a2*x",'Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Spring",'Range','A1');
writematrix("a1",'Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Spring",'Range','A2');
writematrix("a2",'Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Spring",'Range','A3');
writetable(ff.Coefficients,'Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Spring",'Range','B:E');
writematrix('R^2','Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Spring",'Range','F1');
writematrix(ff.Rsquared.Ordinary,'Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Spring",'Range','F2');
elseif Rf==f1.Rsquared.Ordinary
    subplot(2,1,1)
    
errorbar(L(:,i),K(i,:),SDK(i,:),'o');
hold on;
    if Rf<0.5
plot(L(:,i),predict(f1,L(:,i)),L(:,i),mean(double(K(i,:)),'omitnan').*ones(1,length(PL)),'--');
legend(lk,"Fit "+f1.Coefficients{1,1}+" + "+f1.Coefficients{2,1}+" *x", "Mean Spring value "+mean(double(K(i,:)),'omitnan'),'Location','southeast');
    else
        plot(L(:,i),predict(f1,L(:,i)));
        legend(lk,"Fit "+f1.Coefficients{1,1}+" + "+f1.Coefficients{2,1}+" *x",'Location','southeast');
    end
ylabel("K (Nm/rad)");
xlabel("Segment length (mm)");
title(ck);


hold off;
writematrix("y=a1+a2*x",'Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Spring",'Range','A1');
writematrix("a1",'Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Spring",'Range','A2');
writematrix("a2",'Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Spring",'Range','A3');
writetable(f1.Coefficients,'Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Spring",'Range','B:E');
writematrix('R^2','Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Spring",'Range','F1');
writematrix(f1.Rsquared.Ordinary,'Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Spring",'Range','F2');
elseif Rf==f2.Rsquared.Ordinary
 
    subplot(2,1,1)
errorbar(PL,K(i,:),SDK(i,:),'o');
 hold on;
 if Rf<0.5
plot(PL,predict(f2,PL),PL,mean(double(K(i,:)),'omitnan').*ones(1,length(PL)),'--');
legend(lk,"Fit "+f2.Coefficients{1,1}+" + "+f2.Coefficients{2,1}+" *x","Mean Spring value "+mean(double(K(i,:)),'omitnan'),'Location','southeast');
 else
     plot(PL,predict(f2,PL));
     legend(lk,"Fit "+f2.Coefficients{1,1}+" + "+f2.Coefficients{2,1}+" *x",'Location','southeast');
 end
 
ylabel("K (Nm/rad)");
xlabel("Palm length (mm)");
title(ck);


hold off; 
writematrix("y=a1+a2*x",'Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Spring",'Range','A1');
writematrix("a1",'Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Spring",'Range','A2');
writematrix("a2",'Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Spring",'Range','A3');
writetable(f2.Coefficients,'Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Spring",'Range','B:E');
writematrix('R^2','Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Spring",'Range','F1');
writematrix(f2.Rsquared.Ordinary,'Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Spring",'Range','F2');

elseif Rf==f3.Rsquared.Ordinary
  subplot(2,1,1)
 
errorbar(PL.*HB,K(i,:),SDK(i,:),'o');
 hold on;
  if Rf<0.5
plot(PL.*HB,predict(f3,PL.*HB),PL.*HB,mean(double(K(i,:)),'omitnan').*ones(1,length(PL)),'--');
legend(lk,"Fit "+f3.Coefficients{1,1}+" + "+f3.Coefficients{2,1}+" *x","Mean Spring value "+mean(double(K(i,:)),'omitnan'),'Location','southeast');
 else
     plot(PL.*HB,predict(f3,PL.*HB));
     legend(lk,"Fit "+f3.Coefficients{1,1}+" + "+f3.Coefficients{2,1}+" *x",'Location','southeast');
 end
ylabel("K (Nm/rad)");
xlabel("Product of Palm length and Hand breadth (mm^2)");
title(ck);


hold off;  
writematrix("y=a1+a2*x",'Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Spring",'Range','A1');
writematrix("a1",'Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Spring",'Range','A2');
writematrix("a2",'Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Spring",'Range','A3');
writetable(f3.Coefficients,'Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Spring",'Range','B:E');
writematrix('R^2','Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Spring",'Range','F1');
writematrix(f3.Rsquared.Ordinary,'Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Spring",'Range','F2');

elseif Rf==f4.Rsquared.Ordinary 
    
subplot(2,1,1)

errorbar(HLHB,K(i,:),SDK(i,:),'o');
hold on;
if Rf<0.5
plot(HLHB,predict(f4,HLHB),HLHB,mean(double(K(i,:)),'omitnan').*ones(1,length(HLHB)),'--');
legend(lk,"Fit "+f4.Coefficients{1,1}+" + "+f4.Coefficients{2,1}+" *x","Mean Spring value "+mean(double(K(i,:)),'omitnan'),'Location','southeast');
 else
     plot(HLHB,predict(f4,HLHB));
     legend(lk,"Fit "+f4.Coefficients{1,1}+" + "+f4.Coefficients{2,1}+" *x",'Location','southeast');
 end
ylabel("K (Nm/rad)");
xlabel("Product of hand length and Hand breadth (mm^2)");
title(ck);


hold off;  
writematrix("y=a1+a2*x",'Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Spring",'Range','A1');
writematrix("a1",'Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Spring",'Range','A2');
writematrix("a2",'Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Spring",'Range','A3');
writetable(f4.Coefficients,'Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Spring",'Range','B:E');
writematrix('R^2','Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Spring",'Range','F1');
writematrix(f4.Rsquared.Ordinary,'Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Spring",'Range','F2');

elseif Rf==f5.Rsquared.Ordinary 
    
subplot(2,1,1)

errorbar(Ra(:,i),K(i,:),SDK(i,:),'o');
hold on;
if Rf<0.5
plot(Ra(:,i),predict(f5,Ra(:,i)),Ra(:,i),mean(double(K(i,:)),'omitnan').*ones(1,length(HLHB)),'--');
legend(lk,"Fit "+f5.Coefficients{1,1}+" + "+f5.Coefficients{2,1}+" *x","Mean Spring value "+mean(double(K(i,:)),'omitnan'),'Location','southeast');
 else
     plot(Ra(:,i),predict(f5,Ra(:,i)));
     legend(lk,"Fit "+f5.Coefficients{1,1}+" + "+f5.Coefficients{2,1}+" *x",'Location','southeast');
 end
ylabel("K (Nm/rad)");
xlabel("Segment radius (mm)");
title(ck);


hold off;  
writematrix("y=a1+a2*x",'Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Spring",'Range','A1');
writematrix("a1",'Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Spring",'Range','A2');
writematrix("a2",'Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Spring",'Range','A3');
writetable(f5.Coefficients,'Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Spring",'Range','B:E');
writematrix('R^2','Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Spring",'Range','F1');
writematrix(f5.Rsquared.Ordinary,'Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Spring",'Range','F2');

else
subplot(2,1,1)

errorbar(Palm,K(i,:),SDK(i,:),'o');
hold on;
if Rf<0.5
plot(Palm,predict(f6,Palm),Palm,mean(double(K(i,:)),'omitnan').*ones(1,length(HLHB)),'--');
legend(lk,"Fit "+f6.Coefficients{1,1}+" + "+f6.Coefficients{2,1}+" *x","Mean Spring value "+mean(double(K(i,:)),'omitnan'),'Location','southeast');
 else
     plot(Palm,predict(f6,Palm));
     legend(lk,"Fit "+f6.Coefficients{1,1}+" + "+f6.Coefficients{2,1}+" *x",'Location','southeast');
 end
ylabel("K (Nm/rad)");
xlabel("Sum of palm dimensions (mm)");
title(ck);


hold off;  
writematrix("y=a1+a2*x",'Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Spring",'Range','A1');
writematrix("a1",'Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Spring",'Range','A2');
writematrix("a2",'Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Spring",'Range','A3');
writetable(f6.Coefficients,'Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Spring",'Range','B:E');
writematrix('R^2','Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Spring",'Range','F1');
writematrix(f6.Rsquared.Ordinary,'Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Spring",'Range','F2');
       
end

if Rf<0.5
    A=[mean(double(K(i,:)),'omitnan'); std(double(K(i,:)),'omitnan')];
    writematrix("Mean Spring value (Nm/rad)",'Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Spring",'Range','G1');
    writematrix("STD Spring value (Nm/rad)",'Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Spring",'Range','G2');
    writematrix(A,'Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Spring",'Range','H1:H2');
    writematrix(transpose(double(K(i,:))),'Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Spring",'Range','I2');
    writematrix("Spring values",'Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Spring",'Range','I1');
end


if Rff==ff1.Rsquared.Ordinary
subplot(2,1,2)

errorbar(L(:,i).*Ra(:,i),B(i,:),SDB(i,:),'o');
hold on;
if Rff<0.5
plot(L(:,i).*Ra(:,i),predict(ff1,L(:,i).*Ra(:,i)),L(:,i).*Ra(:,i),mean(double(B(i,:)),'omitnan').*ones(1,length(PL)),'--');
legend(ld,"Fit "+ff1.Coefficients{1,1}+" + "+ff1.Coefficients{2,1}+" *x","Mean Damper value "+mean(double(B(i,:)),'omitnan'),'Location','southeast');
 else
     plot(L(:,i).*Ra(:,i),predict(ff1,L(:,i).*Ra(:,i)));
     legend(ld,"Fit "+ff1.Coefficients{1,1}+" + "+ff1.Coefficients{2,1}+" *x",'Location','southeast');
 end
ylabel("B (Nms/rad)");
xlabel("Product of segment length and radius (mm^2)");
title(cb);


hold off;
writematrix("y=a1+a2*x",'Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Damper",'Range','A1');
writematrix("a1",'Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Damper",'Range','A2');
writematrix("a2",'Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Damper",'Range','A3');
writetable(ff1.Coefficients,'Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Damper",'Range','B:E');
writematrix('R^2','Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Damper",'Range','F1');
writematrix(ff1.Rsquared.Ordinary,'Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Damper",'Range','F2');

elseif Rff==ff2.Rsquared.Ordinary
    subplot(2,1,2)
   
errorbar(Ra(:,i),B(i,:),SDB(i,:),'o');
 hold on;
  if Rff<0.5
plot(Ra(:,i),predict(ff2,Ra(:,i)),Ra(:,i),mean(double(B(i,:)),'omitnan').*ones(1,length(PL)),'--');
legend(ld,"Fit "+ff2.Coefficients{1,1}+" + "+ff2.Coefficients{2,1}+" *x","Mean Damper value "+mean(double(B(i,:)),'omitnan'),'Location','southeast');
 else
     plot(Ra(:,i),predict(ff2,Ra(:,i)));
     legend(ld,"Fit "+ff2.Coefficients{1,1}+" + "+ff2.Coefficients{2,1}+" *x",'Location','southeast');
 end
ylabel("B (Nms/rad)");
xlabel("Segment radius (mm)");
title(cb);


hold off;
writematrix("y=a1+a2*x",'Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Damper",'Range','A1');
writematrix("a1",'Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Damper",'Range','A2');
writematrix("a2",'Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Damper",'Range','A3');
writetable(ff2.Coefficients,'Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Damper",'Range','B:E');
writematrix('R^2','Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Damper",'Range','F1');
writematrix(ff2.Rsquared.Ordinary,'Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Damper",'Range','F2');

elseif Rff==ff3.Rsquared.Ordinary
 subplot(2,1,2)

errorbar(PL,B(i,:),SDB(i,:),'o');
 hold on;
 if Rff<0.5
plot(PL,predict(ff3,PL),PL,mean(double(B(i,:)),'omitnan').*ones(1,length(PL)),'--');
legend(ld,"Fit "+ff3.Coefficients{1,1}+" + "+ff3.Coefficients{2,1}+" *x","Mean Damper value "+mean(double(B(i,:)),'omitnan'),'Location','southeast');
 else
     plot(PL,predict(ff3,PL));
     legend(ld,"Fit "+ff3.Coefficients{1,1}+" + "+ff3.Coefficients{2,1}+" *x",'Location','southeast');
 end
ylabel("B (Nms/rad)");
xlabel("Palm length (mm)");
title(cb);


hold off; 
writematrix("y=a1+a2*x",'Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Damper",'Range','A1');
writematrix("a1",'Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Damper",'Range','A2');
writematrix("a2",'Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Damper",'Range','A3');
writetable(ff3.Coefficients,'Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Damper",'Range','B:E');
writematrix('R^2','Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Damper",'Range','F1');
writematrix(ff3.Rsquared.Ordinary,'Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Damper",'Range','F2');

elseif Rff==ff4.Rsquared.Ordinary
  subplot(2,1,2)
  
errorbar(PL.*HB,B(i,:),SDB(i,:),'o');
hold on;
  if Rff<0.5
plot(PL.*HB,predict(ff4,PL.*HB),PL.*HB,mean(double(B(i,:)),'omitnan').*ones(1,length(PL)),'--');
legend(ld,"Fit "+ff4.Coefficients{1,1}+" + "+ff4.Coefficients{2,1}+" *x","Mean Damper value "+mean(double(B(i,:)),'omitnan'),'Location','southeast');
 else
     plot(PL.*HB,predict(ff4,PL.*HB));
     legend(ld,"Fit "+ff4.Coefficients{1,1}+" + "+ff4.Coefficients{2,1}+" *x",'Location','southeast');
 end
ylabel("B (Nms/rad)");
xlabel("Product of Palm length and Hand breadth (mm^2)");
title(cb);


hold off;  
writematrix("y=a1+a2*x",'Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Damper",'Range','A1');
writematrix("a1",'Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Damper",'Range','A2');
writematrix("a2",'Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Damper",'Range','A3');
writetable(ff4.Coefficients,'Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Damper",'Range','B:E');
writematrix('R^2','Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Damper",'Range','F1');
writematrix(ff4.Rsquared.Ordinary,'Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Damper",'Range','F2');

elseif Rff==ff5.Rsquared.Ordinary 
    
subplot(2,1,2)

errorbar(HLHB,B(i,:),SDB(i,:),'o');
hold on;
if Rff<0.5
plot(HLHB,predict(ff5,HLHB),HLHB,mean(double(B(i,:)),'omitnan').*ones(1,length(HLHB)),'--');
legend(ld,"Fit "+ff5.Coefficients{1,1}+" + "+ff5.Coefficients{2,1}+" *x","Mean Damper value "+mean(double(B(i,:)),'omitnan'),'Location','southeast');
 else
     plot(HLHB,predict(ff5,HLHB));
     legend(ld,"Fit "+ff5.Coefficients{1,1}+" + "+ff5.Coefficients{2,1}+" *x",'Location','southeast');
 end
ylabel("B (Nms/rad)");
xlabel("Product of hand length and Hand breadth (mm^2)");
title(cb);


hold off;  
writematrix("y=a1+a2*x",'Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Damper",'Range','A1');
writematrix("a1",'Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Damper",'Range','A2');
writematrix("a2",'Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Damper",'Range','A3');
writetable(ff5.Coefficients,'Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Damper",'Range','B:E');
writematrix('R^2','Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Damper",'Range','F1');
writematrix(ff5.Rsquared.Ordinary,'Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Damper",'Range','F2');

elseif Rff==ff6.Rsquared.Ordinary

    subplot(2,1,2)

errorbar(L(:,i),B(i,:),SDB(i,:),'o');
hold on;
if Rff<0.5
plot(L(:,i),predict(ff6,L(:,i)),L(:,i),mean(double(B(i,:)),'omitnan').*ones(1,length(HLHB)),'--');
legend(ld,"Fit "+ff6.Coefficients{1,1}+" + "+ff6.Coefficients{2,1}+" *x","Mean Damper value "+mean(double(B(i,:)),'omitnan'),'Location','southeast');
 else
     plot(L(:,i),predict(ff6,L(:,i)));
     legend(ld,"Fit "+ff6.Coefficients{1,1}+" + "+ff6.Coefficients{2,1}+" *x",'Location','southeast');
 end
ylabel("B (Nms/rad)");
xlabel("Segment length (mm)");
title(cb);


hold off;  
writematrix("y=a1+a2*x",'Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Damper",'Range','A1');
writematrix("a1",'Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Damper",'Range','A2');
writematrix("a2",'Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Damper",'Range','A3');
writetable(ff6.Coefficients,'Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Damper",'Range','B:E');
writematrix('R^2','Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Damper",'Range','F1');
writematrix(ff6.Rsquared.Ordinary,'Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Damper",'Range','F2');

else
    subplot(2,1,2)

errorbar(Palm,B(i,:),SDB(i,:),'o');
hold on;
if Rff<0.5
plot(Palm,predict(ff7,Palm),Palm,mean(double(B(i,:)),'omitnan').*ones(1,length(HLHB)),'--');
legend(ld,"Fit "+ff7.Coefficients{1,1}+" + "+ff7.Coefficients{2,1}+" *x","Mean Damper value "+mean(double(B(i,:)),'omitnan'),'Location','southeast');
 else
     plot(Palm,predict(ff7,Palm));
     legend(ld,"Fit "+ff7.Coefficients{1,1}+" + "+ff7.Coefficients{2,1}+" *x",'Location','southeast');
 end
ylabel("B (Nms/rad)");
xlabel("Sum of palm dimensions (mm)");
title(cb);


hold off;  
writematrix("y=a1+a2*x",'Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Damper",'Range','A1');
writematrix("a1",'Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Damper",'Range','A2');
writematrix("a2",'Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Damper",'Range','A3');
writetable(ff7.Coefficients,'Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Damper",'Range','B:E');
writematrix('R^2','Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Damper",'Range','F1');
writematrix(ff7.Rsquared.Ordinary,'Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Damper",'Range','F2');

end
if Rff<0.5
     writematrix("Mean Damper value (Nms/rad)",'Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Damper",'Range','G1');
    writematrix("STD Damper value (Nms/rad)",'Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Damper",'Range','G2');
    writematrix([mean(double(B(i,:)),'omitnan'); std(double(B(i,:)),'omitnan')],'Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Damper",'Range','H1:H2');
    writematrix(transpose(double(B(i,:))),'Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Damper",'Range','I2');
    writematrix("Damper values",'Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Damper",'Range','I1');
end
saveas(h,cp,'jpg');
end
Ls=[Lt(:,1)+Lt(:,2)+Lt(:,3) Li(:,1)+Li(:,2)+Li(:,3) Lm(:,1)+Lm(:,2)+Lm(:,3) Lr(:,1)+Lr(:,2)+Lr(:,3) Ll(:,1)+Ll(:,2)+Ll(:,3)];
Rs=[Rt(:,1)+Rt(:,2)+Rt(:,3) Ri(:,1)+Ri(:,2)+Ri(:,3) Rm(:,1)+Rm(:,2)+Rm(:,3) Rr(:,1)+Rr(:,2)+Rr(:,3) Rl(:,1)+Rl(:,2)+Rl(:,3)];

%Abduction movment scaling functions
for s=1:5
w=1./SDKa(s,:);
wb=1./SDBa(s,:);
wb=double(wb./max(wb));
wk=double(w./max(w));

%Spring scaling function
ff=fitnlm(La(:,s),double(Ka(s,:)),model,[0.003 0.01],'Weight',transpose(wk),'Options',opts);
f2=fitnlm(Ls(:,s),double(Ka(s,:)),model,[0.003 0.01],'Weight',transpose(wk),'Options',opts);
f3=fitnlm(PL,double(Ka(s,:)),model,[0.003 0.01],'Weight',transpose(wk),'Options',opts);
f4=fitnlm(PLHB,double(Ka(s,:)),model,[0.003 0.01],'Weight',transpose(wk),'Options',opts);
f5=fitnlm(HLHB,double(Ka(s,:)),model,[0.003 0.01],'Weight',transpose(wk),'Options',opts);
f6=fitnlm(Rs(:,s),double(Ka(s,:)),model,[0.003 0.01],'Weight',transpose(wk),'Options',opts);
f7=fitnlm(Palm,double(Ka(s,:)),model,[0.003 0.01],'Weight',transpose(wk),'Options',opts);

%Damper scaling functions
ff1=fitnlm(La(:,s),double(Ba(s,:)),model,[0.003 0.01],'Weight',transpose(wb),'Options',opts);
ff2=fitnlm(Rs(:,s),double(Ba(s,:)),model,[0.003 0.01],'Weight',transpose(wb),'Options',opts);
ff3=fitnlm(PL,double(Ba(s,:)),model,[0.003 0.01],'Weight',transpose(wb),'Options',opts);
ff4=fitnlm(PLHB,double(Ba(s,:)),model,[0.003 0.01],'Weight',transpose(wb),'Options',opts);
ff5=fitnlm(HLHB,double(Ba(s,:)),model,[0.003 0.01],'Weight',transpose(wb),'Options',opts);
ff6=fitnlm(Ls(:,s),double(Ba(s,:)),model,[0.003 0.01],'Weight',transpose(wb),'Options',opts);
ff7=fitnlm(Palm,double(Ba(s,:)),model,[0.003 0.01],'Weight',transpose(wb),'Options',opts);

%The R^2 parameters from NLM
[Rf,kl]=max([ff.Rsquared.Ordinary  f2.Rsquared.Ordinary f3.Rsquared.Ordinary f4.Rsquared.Ordinary f5.Rsquared.Ordinary f6.Rsquared.Ordinary f7.Rsquared.Ordinary]);
[Rff,kl1]=max([ff1.Rsquared.Ordinary ff2.Rsquared.Ordinary ff3.Rsquared.Ordinary ff4.Rsquared.Ordinary ff5.Rsquared.Ordinary ff6.Rsquared.Ordinary ff7.Rsquared.Ordinary ]);

h=figure('Position',get(0,'Screensize'));
if s==1
    p="thumb";
elseif s==2
    p="index";
elseif s==3
    p="middle";
    elseif s==4
        p="ring";
else
    p="little";
end
ShName=p+" Abduction";

if Rf== ff.Rsquared.Ordinary 
subplot(2,1,1)

errorbar(La(:,s),Ka(s,:),SDKa(s,:),'o');
hold on;
if Rf<0.5
plot(La(:,s),predict(ff,La(:,s)),La(:,s),mean(double(Ka(s,:)),'omitnan').*ones(1,length(PL)),'--');
legend("Spring Abduction data of "+p,"Fit "+ff.Coefficients{1,1}+" + "+ff.Coefficients{2,1}+" *x","Mean Spring value "+mean(double(Ka(s,:)),'omitnan'),'Location','southeast');
 else
     plot(La(:,s),predict(ff,La(:,s)));
     legend("Spring Abduction data of "+p,"Fit "+ff.Coefficients{1,1}+" + "+ff.Coefficients{2,1}+" *x",'Location','southeast');
end
ylabel("K (Nm/rad)");
xlabel("Sum of Products of segment lengths and radii (mm^2)");
title("Abduction Spring constant regression  of "+p+" finger with R^2 "+Rf+" and RMSE "+ff.RMSE+" Nm/rad");


hold off;
writematrix("y=a1+a2*x",'Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Spring",'Range','A1');
writematrix("a1",'Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Spring",'Range','A2');
writematrix("a2",'Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Spring",'Range','A3');
writetable(ff.Coefficients,'Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Spring");
writematrix('R^2','Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Spring",'Range','F1');
writematrix(ff.Rsquared.Ordinary,'Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Spring",'Range','F2');

elseif Rf==f2.Rsquared.Ordinary
  subplot(2,1,1)

errorbar(Ls(:,s),Ka(s,:),SDKa(s,:),'o');
  hold on;
  if Rf<0.5
plot(Ls(:,s),predict(f2,Ls(:,s)),Ls(:,s),mean(double(Ka(s,:)),'omitnan').*ones(1,length(PL)),'--');
legend("Spring Abduction data of "+p,"Fit "+f2.Coefficients{1,1}+" + "+f2.Coefficients{2,1}+" *x","Mean Spring value "+mean(double(Ka(s,:)),'omitnan'),'Location','southeast');
 else
     plot(Ls(:,s),predict(f2,Ls(:,s)));
     legend("Spring Abduction data of "+p,"Fit "+f2.Coefficients{1,1}+" + "+f2.Coefficients{2,1}+" *x",'Location','southeast');
 end
ylabel("K (Nm/rad)");
xlabel("Sum of segment lengths (mm)");
title("Abduction Spring constant regression  of "+p+" finger with R^2 "+Rf+" and RMSE "+f2.RMSE+" Nm/rad");


hold off;  
writematrix("y=a1+a2*x",'Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Spring",'Range','A1');
writematrix("a1",'Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Spring",'Range','A2');
writematrix("a2",'Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Spring",'Range','A3');
writetable(f2.Coefficients,'Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Spring");
writematrix('R^2','Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Spring",'Range','F1');
writematrix(f2.Rsquared.Ordinary,'Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Spring",'Range','F2');

elseif Rf==f3.Rsquared.Ordinary
subplot(2,1,1)

errorbar(PL,Ka(s,:),SDKa(s,:),'o');
hold on;
if Rf<0.5
plot(PL,predict(f3,PL),PL,mean(double(Ka(s,:)),'omitnan').*ones(1,length(PL)),'--');
legend("Spring Abduction data of "+p,"Fit "+f3.Coefficients{1,1}+" + "+f3.Coefficients{2,1}+" *x","Mean Spring value "+mean(double(Ka(s,:)),'omitnan'),'Location','southeast');
 else
     plot(PL,predict(f3,PL));
     legend("Spring Abduction data of "+p,"Fit "+f3.Coefficients{1,1}+" + "+f3.Coefficients{2,1}+" *x",'Location','southeast');
 end
ylabel("K (Nm/rad)");
xlabel("Palm Length (mm)");
title("Abduction Spring constant regression  of "+p+" finger with R^2 "+Rf+" and RMSE "+f3.RMSE+" Nm/rad");


hold off;   
writematrix("y=a1+a2*x",'Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Spring",'Range','A1');
writematrix("a1",'Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Spring",'Range','A2');
writematrix("a2",'Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Spring",'Range','A3');
writetable(f3.Coefficients,'Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Spring");
writematrix('R^2','Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Spring",'Range','F1');
writematrix(f3.Rsquared.Ordinary,'Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Spring",'Range','F2');

elseif Rf==f4.Rsquared.Ordinary
    subplot(2,1,1)
  
errorbar(PLHB,Ka(s,:),SDKa(s,:),'o');
  hold on;
 if Rf<0.5
plot(PLHB,predict(f4,PLHB),PLHB,mean(double(Ka(s,:)),'omitnan').*ones(1,length(PLHB)),'--');
legend("Spring Abduction data of "+p,"Fit "+f4.Coefficients{1,1}+" + "+f4.Coefficients{2,1}+" *x","Mean Spring value "+mean(double(Ka(s,:)),'omitnan'),'Location','southeast');
 else
     plot(PLHB,predict(f4,PLHB));
     legend("Spring Abduction data of "+p,"Fit "+f4.Coefficients{1,1}+" + "+f4.Coefficients{2,1}+" *x",'Location','southeast');
 end
ylabel("K (Nm/rad)");
xlabel("Product of Palm length and Hand breadth (mm^2)");
title("Abduction Spring constant regression  of "+p+" finger with R^2 "+Rf+ "and RMSE "+f4.RMSE+" Nm/rad");


hold off;  
writematrix("y=a1+a2*x",'Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Spring",'Range','A1');
writematrix("a1",'Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Spring",'Range','A2');
writematrix("a2",'Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Spring",'Range','A3');
writetable(f4.Coefficients,'Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Spring");
writematrix('R^2','Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Spring",'Range','F1');
writematrix(f4.Rsquared.Ordinary,'Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Spring",'Range','F2');

elseif Rf==f5.Rsquared.Ordinary
  subplot(2,1,1)
 
errorbar(HLHB,Ka(s,:),SDKa(s,:),'o');
 hold on;
  if Rf<0.5
plot(HLHB,predict(f5,HLHB),HLHB,mean(double(Ka(s,:)),'omitnan').*ones(1,length(HLHB)),'--');
legend("Spring Abduction data of "+p,"Fit "+f5.Coefficients{1,1}+" + "+f5.Coefficients{2,1}+" *x","Mean Spring value "+mean(double(Ka(s,:)),'omitnan'),'Location','southeast');
 else
     plot(HLHB,predict(f5,HLHB));
     legend("Spring Abduction data of "+p,"Fit "+f5.Coefficients{1,1}+" + "+f5.Coefficients{2,1}+" *x",'Location','southeast');
 end
ylabel("K (Nm/rad)");
xlabel("Product of Hand length and Hand breadth (mm^2)");
title("Abduction Spring constant regression  of "+p+" finger with R^2 "+Rf+" and RMSE "+f5.RMSE+" Nm/rad");


hold off;  
writematrix("y=a1+a2*x",'Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Spring",'Range','A1');
writematrix("a1",'Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Spring",'Range','A2');
writematrix("a2",'Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Spring",'Range','A3');
writetable(f5.Coefficients,'Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Spring");
writematrix('R^2','Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Spring",'Range','F1');
writematrix(f5.Rsquared.Ordinary,'Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Spring",'Range','F2');

elseif Rf==f6.Rsquared.Ordinary
    subplot(2,1,1)
 
errorbar(Rs(:,s),Ka(s,:),SDKa(s,:),'o');
 hold on;
  if Rf<0.5
plot(Rs(:,s),predict(f6,Rs(:,s)),Rs(:,s),mean(double(Ka(s,:)),'omitnan').*ones(1,length(HLHB)),'--');
legend("Spring Abduction data of "+p,"Fit "+f6.Coefficients{1,1}+" + "+f6.Coefficients{2,1}+" *x","Mean Spring value "+mean(double(Ka(s,:)),'omitnan'),'Location','southeast');
 else
     plot(Rs(:,s),predict(f6,Rs(:,s)));
     legend("Spring Abduction data of "+p,"Fit "+f6.Coefficients{1,1}+" + "+f6.Coefficients{2,1}+" *x",'Location','southeast');
 end
ylabel("K (Nm/rad)");
xlabel("Sum of segment radii (mm)");
title("Abduction Spring constant regression  of "+p+" finger with R^2 "+Rf+" and RMSE "+f6.RMSE+" Nm/rad");

hold off;  
writematrix("y=a1+a2*x",'Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Spring",'Range','A1');
writematrix("a1",'Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Spring",'Range','A2');
writematrix("a2",'Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Spring",'Range','A3');
writetable(f6.Coefficients,'Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Spring");
writematrix('R^2','Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Spring",'Range','F1');
writematrix(f6.Rsquared.Ordinary,'Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Spring",'Range','F2');

else
    subplot(2,1,1)
 
errorbar(Palm,Ka(s,:),SDKa(s,:),'o');
 hold on;
  if Rf<0.5
plot(Palm,predict(f7,Palm),Palm,mean(double(Ka(s,:)),'omitnan').*ones(1,length(HLHB)),'--');
legend("Spring Abduction data of "+p,"Fit "+f7.Coefficients{1,1}+" + "+f7.Coefficients{2,1}+" *x","Mean Spring value "+mean(double(Ka(s,:)),'omitnan'),'Location','southeast');
 else
     plot(Palm,predict(f7,Palm));
     legend("Spring Abduction data of "+p,"Fit "+f7.Coefficients{1,1}+" + "+f7.Coefficients{2,1}+" *x",'Location','southeast');
 end
ylabel("K (Nm/rad)");
xlabel("Sum of palm dimensions (mm)");
title("Abduction Spring constant regression  of "+p+" finger with R^2 "+Rf+" and RMSE "+f7.RMSE+" Nm/rad");


hold off;  
writematrix("y=a1+a2*x",'Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Spring",'Range','A1');
writematrix("a1",'Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Spring",'Range','A2');
writematrix("a2",'Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Spring",'Range','A3');
writetable(f7.Coefficients,'Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Spring");
writematrix('R^2','Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Spring",'Range','F1');
writematrix(f7.Rsquared.Ordinary,'Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Spring",'Range','F2');
end

if Rf<0.5
  writematrix("Mean Spring value (Nm/rad)",'Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Spring",'Range','G1');
    writematrix("STD Spring value (Nm/rad)",'Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Spring",'Range','G2');
    writematrix([mean(double(Ka(s,:)),'omitnan'); std(double(Ka(s,:)),'omitnan')],'Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Spring",'Range','H1:H2'); 
    writematrix(transpose(double(Ka(s,:))),'Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Spring",'Range','I2');
    writematrix("Spring values",'Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Spring",'Range','I1');
end


if Rff==ff1.Rsquared.Ordinary
subplot(2,1,2)

errorbar(La(:,s),Ba(s,:),SDBa(s,:),'o');
hold on;
if Rff<0.5
plot(La(:,s),predict(ff1,La(:,s)),La(:,s),mean(double(Ba(s,:)),'omitnan').*ones(1,length(PL)),'--');
legend("Damper Abduction data of "+p,"Fit "+ff1.Coefficients{1,1}+" + "+ff1.Coefficients{2,1}+" *x","Mean Damper value "+mean(double(Ba(s,:)),'omitnan'),'Location','southeast');
 else
     plot(La(:,s),predict(ff1,La(:,s)));
     legend("Damper Abduction data of "+p,"Fit "+ff1.Coefficients{1,1}+" + "+ff1.Coefficients{2,1}+" *x",'Location','southeast');
 end
ylabel("B (Nms/rad)");
xlabel("Sum of Product of segment lengths and radii (mm^2)");
title("Abduction Damper constant regression  of "+p+" finger with R^2 "+Rff+" and RMSE "+ff1.RMSE+" Nms/rad");


hold off;
writematrix("y=a1+a2*x",'Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Damper",'Range','A1');
writematrix("a1",'Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Damper",'Range','A2');
writematrix("a2",'Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Damper",'Range','A3');
writetable(ff1.Coefficients,'Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Damper");
writematrix('R^2','Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Damper",'Range','F1');
writematrix(ff1.Rsquared.Ordinary,'Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Damper",'Range','F2');

elseif Rff==ff2.Rsquared.Ordinary
 subplot(2,1,2)
 
errorbar(Rs(:,s),Ba(s,:),SDBa(s,:),'o');
hold on;
 if Rff<0.5
plot(Rs(:,s),predict(ff2,Rs(:,s)),Rs(:,s),mean(Ba(s,:)).*ones(1,length(PL)),'--');
legend("Damper Abduction data of "+p,"Fit "+ff2.Coefficients{1,1}+" + "+ff2.Coefficients{2,1}+" *x","Mean Damper value "+mean(double(Ba(s,:)),'omitnan'),'Location','southeast');
 else
     plot(Rs(:,s),predict(ff2,Rs(:,s)));
     legend("Damper Abduction data of "+p,"Fit "+ff2.Coefficients{1,1}+" + "+ff2.Coefficients{2,1}+" *x",'Location','southeast');
 end
ylabel("B (Nms/rad)");
xlabel("Sum of segment radii (mm)");
title("Abduction Damper constant regression  of "+p+" finger with R^2 "+Rff+" and RMSE "+ff2.RMSE+" Nms/rad");


hold off;   
writematrix("y=a1+a2*x",'Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Damper",'Range','A1');
writematrix("a1",'Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Damper",'Range','A2');
writematrix("a2",'Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Damper",'Range','A3');
writetable(ff2.Coefficients,'Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Damper");
writematrix('R^2','Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Damper",'Range','F1');
writematrix(ff2.Rsquared.Ordinary,'Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Damper",'Range','F2');

elseif Rff==ff3.Rsquared.Ordinary
    subplot(2,1,2)

errorbar(PL,Ba(s,:),SDBa(s,:),'o');
    hold on;
 if Rff<0.5
plot(PL,predict(ff3,PL),PL,mean(Ba(s,:)).*ones(1,length(Palm)),'--');
legend("Damper Abduction data of "+p,"Fit "+ff3.Coefficients{1,1}+" + "+ff3.Coefficients{2,1}+" *x","Mean Damper value "+mean(double(Ba(s,:)),'omitnan'),'Location','southeast');
 else
     plot(PL,predict(ff3,PL));
     legend("Damper Abduction data of "+p,"Fit "+ff3.Coefficients{1,1}+" + "+ff3.Coefficients{2,1}+" *x",'Location','southeast');
 end
ylabel("B (Nms/rad)");
xlabel("Palm Length (mm)");
title("Abduction Damper constant regression  of "+p+" finger with R^2 "+Rff+" and RMSE "+ff3.RMSE+" Nms/rad");


hold off;   
writematrix("y=a1+a2*x",'Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Damper",'Range','A1');
writematrix("a1",'Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Damper",'Range','A2');
writematrix("a2",'Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Damper",'Range','A3');
writetable(ff3.Coefficients,'Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Damper");  
writematrix('R^2','Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Damper",'Range','F1');
writematrix(ff3.Rsquared.Ordinary,'Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Damper",'Range','F2');

elseif Rff==ff4.Rsquared.Ordinary
    subplot(2,1,2)
  
errorbar(PLHB,Ba(s,:),SDBa(s,:),'o');
  hold on;
 if Rff<0.5
plot(PLHB,predict(ff4,PLHB),PLHB,mean(Ba(s,:)).*ones(1,length(PLHB)),'--');
legend("Damper Abduction data of "+p,"Fit "+ff4.Coefficients{1,1}+" + "+ff4.Coefficients{2,1}+" *x","Mean Damper value "+mean(double(Ba(s,:)),'omitnan'),'Location','southeast');
 else
     plot(PLHB,predict(ff4,PLHB));
     legend("Damper Abduction data of "+p,"Fit "+ff4.Coefficients{1,1}+" + "+ff4.Coefficients{2,1}+" *x",'Location','southeast');
 end
ylabel("B (Nms/rad)");
xlabel("Product of Palm length and hand breadth (mm^2)");
title("Abduction Damper constant regression  of "+p+" finger with R^2 "+Rff+" and RMSE "+ff4.RMSE+" Nms/rad");


hold off; 
writematrix("y=a1+a2*x",'Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Damper",'Range','A1');
writematrix("a1",'Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Damper",'Range','A2');
writematrix("a2",'Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Damper",'Range','A3');
writetable(ff4.Coefficients,'Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Damper");
writematrix('R^2','Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Damper",'Range','F1');
writematrix(ff4.Rsquared.Ordinary,'Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Damper",'Range','F2');

elseif Rff==ff5.Rsquared.Ordinary
    subplot(2,1,2)
    
errorbar(HLHB,Ba(s,:),SDBa(s,:),'o');
hold on;
 if Rff<0.5
plot(HLHB,predict(ff5,HLHB),HLHB,mean(Ba(s,:)).*ones(1,length(HLHB)),'--');
legend("Damper Abduction data of "+p,"Fit "+ff5.Coefficients{1,1}+" + "+ff5.Coefficients{2,1}+" *x","Mean Damper value "+mean(double(Ba(s,:)),'omitnan'),'Location','southeast');
 else
     plot(HLHB,predict(ff5,HLHB));
     legend("Damper Abduction data of "+p,"Fit "+ff5.Coefficients{1,1}+" + "+ff5.Coefficients{2,1}+" *x",'Location','southeast');
 end
ylabel("B (Nms/rad)");
xlabel("Product of Hand length and Hand breadth (mm^2)");
title("Abduction Damper constant regression  of "+p+" finger with R^2 "+Rff+" and RMSE "+ff5.RMSE+" Nms/rad");

hold off; 
writematrix("y=a1+a2*x",'Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Damper",'Range','A1');
writematrix("a1",'Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Damper",'Range','A2');
writematrix("a2",'Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Damper",'Range','A3');
writetable(ff5.Coefficients,'Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Damper");
writematrix('R^2','Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Damper",'Range','F1');
writematrix(ff5.Rsquared.Ordinary,'Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Damper",'Range','F2');

elseif Rff==ff6.Rsquared.Ordinary
    subplot(2,1,2)
    
errorbar(Ls(:,s),Ba(s,:),SDBa(s,:),'o');
hold on;
 if Rff<0.5
plot(Ls(:,s),predict(ff6,Ls(:,s)),Ls(:,s),mean(Ba(s,:)).*ones(1,length(HLHB)),'--');
legend("Damper Abduction data of "+p,"Fit "+ff6.Coefficients{1,1}+" + "+ff6.Coefficients{2,1}+" *x","Mean Damper value "+mean(double(Ba(s,:)),'omitnan'),'Location','southeast');
 else
     plot(Ls(:,s),predict(ff6,Ls(:,s)));
     legend("Damper Abduction data of "+p,"Fit "+ff6.Coefficients{1,1}+" + "+ff6.Coefficients{2,1}+" *x",'Location','southeast');
 end
ylabel("B (Nms/rad)");
xlabel("Sum of segment lengths (mm)");
title("Abduction Damper constant regression  of "+p+" finger with R^2 "+Rff+" and RMSE "+ff6.RMSE+" Nms/rad");

hold off; 
writematrix("y=a1+a2*x",'Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Damper",'Range','A1');
writematrix("a1",'Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Damper",'Range','A2');
writematrix("a2",'Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Damper",'Range','A3');
writetable(ff6.Coefficients,'Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Damper");
writematrix('R^2','Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Damper",'Range','F1');
writematrix(ff6.Rsquared.Ordinary,'Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Damper",'Range','F2');

else
    subplot(2,1,2)
    
errorbar(Palm,Ba(s,:),SDBa(s,:),'o');
hold on;
 if Rff<0.5
plot(Palm,predict(ff7,Palm),Palm,mean(double(Ba(s,:)),'omitnan').*ones(1,length(HLHB)),'--');
legend("Damper Abduction data of "+p,"Fit "+ff7.Coefficients{1,1}+" + "+ff7.Coefficients{2,1}+" *x","Mean Damper value "+mean(double(Ba(s,:)),'omitnan'),'Location','southeast');
 else
     plot(Palm,predict(ff7,Palm));
     legend("Damper Abduction data of "+p,"Fit "+ff7.Coefficients{1,1}+" + "+ff7.Coefficients{2,1}+" *x",'Location','southeast');
 end
ylabel("B (Nms/rad)");
xlabel("Sum of palm dimensions (mm)");
title("Abduction Damper constant regression  of "+p+" finger with R^2 "+Rff+" and RMSE "+ff7.RMSE+" Nms/rad");

hold off; 
writematrix("y=a1+a2*x",'Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Damper",'Range','A1');
writematrix("a1",'Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Damper",'Range','A2');
writematrix("a2",'Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Damper",'Range','A3');
writetable(ff7.Coefficients,'Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Damper");
writematrix('R^2','Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Damper",'Range','F1');
writematrix(ff7.Rsquared.Ordinary,'Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Damper",'Range','F2');

end

if Rff<0.5
  writematrix("Mean Damper value (Nms/rad)",'Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Damper",'Range','G1');
    writematrix("STD Damper value (Nms/rad)",'Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Damper",'Range','G2');
    writematrix([mean(double(Ba(s,:)),'omitnan'); std(double(Ba(s,:)),'omitnan')],'Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Damper",'Range','H1:H2');
    writematrix(transpose(double(Ba(s,:))),'Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Damper",'Range','I2');
    writematrix("Damper values",'Non-Linear model fit coefficients.xlsx','Sheet',ShName+" Damper",'Range','I1');
end

saveas(h,"Non linear fit of "+p+" abduction",'jpg');
end