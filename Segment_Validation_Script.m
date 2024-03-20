%This part calculates the difference between estimating the segment length
%using the scaling functions determined here and the ones from Buchholz

%Segment dimensions are in mm and scaling functions take dimensions in mm

clear all;
close all;
clc;

%Set alpha to 0.05 (5%)
a=0.05;

D=["Thumb","Index","Middle","Ring","Little"];
G=input("are you working from your laptop or uni pc?: ",'s');

if G=='l'
    File="C:\Users\panos\OneDrive - University of Warwick\PhD\Hand Trials\Data for Scaling functions\Validation_Set.xlsx";
else
%If i am working from the University PC
File="C:\Users\u1857308\OneDrive - University of Warwick\PhD\Hand Trials\Data for Scaling functions\Validation_Set.xlsx";
end

Thumb=importfile(File,D(1));
Index=importfile(File,D(2));
Middle=importfile(File,D(3));
Ring=importfile(File,D(4));
Little=importfile(File,D(5));
Palm=importfile(File,'Palm');

%Palm Length
PL=Palm.data(:,1);
%Hand Breadth
HB=Palm.data(:,2);
%Hand Width
HW=Palm.data(:,3);

Lt=[Thumb.data(:,5) Thumb.data(:,3) Thumb.data(:,1)];
Li=[Index.data(:,5) Index.data(:,3) Index.data(:,1)];
Lm=[Middle.data(:,5) Middle.data(:,3) Middle.data(:,1)];
Lr=[Ring.data(:,5) Ring.data(:,3) Ring.data(:,1)];
Ll=[Little.data(:,5) Little.data(:,3) Little.data(:,1)];

Rt=[Thumb.data(:,6) Thumb.data(:,4) Thumb.data(:,2)];
Ri=[Index.data(:,6) Index.data(:,4) Index.data(:,2)];
Rm=[Middle.data(:,6) Middle.data(:,4) Middle.data(:,2)];
Rr=[Ring.data(:,6) Ring.data(:,4) Ring.data(:,2)];
Rl=[Little.data(:,6) Little.data(:,4) Little.data(:,2)];

L=[Lt Li Lm Lr Ll];
R=[Rt Ri Rm Rr Rl];

Z=zeros(10,1);
for i=1:10
    s=0;
    for j=1:3
        s=s+Lm(i,j);
    end
    Z(i)=s;
end

%Hand Length
HL=PL+Z;



%Buchholz scaling functions
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


%My scaling functions for segment length
syms y z 

LMt=feval(matlabFunction(25.746+0.5585*z),HW);
LPt=feval(matlabFunction(12.1963+0.3834*x-0.5422*y-0.0073*z),HL,HB,HW);
LDt=feval(matlabFunction(-12.8766+0.2269*x+0.03412*y-0.03416*z),HL,HB,HW);
LPi=feval(matlabFunction(28.3209-0.2144*x+0.3*y+0.087*z),PL,HB,HW);
LMi=feval(matlabFunction(14.693+0.4022*z),HW);
LDi=feval(matlabFunction(-3.2591+0.1106*x+0.0466*y+0.1239*z),HL,HB,HW);
LPm=feval(matlabFunction(12.9316-0.124*x+0.4105*y+0.0677*z),PL,HB,HW);
LMm=feval(matlabFunction(21.6329-0.0767*x+0.1359*y+0.1623*z),PL,HB,HW);
LDm=feval(matlabFunction(8.9697-0.0149*x+0.1297*y+0.2301*z),PL,HB,HW);
LPr=feval(matlabFunction(-11.003+0.3169*x+0.0494*y-0.4805*z),HL,HB,HW);
LMr=feval(matlabFunction(8.0021+0.1227*x+0.0362*y-0.129*z),HL,HB,HW);
LDr=feval(matlabFunction(-6.229+0.1586*x+0.0316*y+0.0018*z),HL,HB,HW);
LPl=feval(matlabFunction(10.264-0.1354*x+0.4077*y-0.1548*z),PL,HB,HW);
LMl=feval(matlabFunction(6.046+0.1015*x-0.018*y-0.09*z),HL,HB,HW);
LDl=feval(matlabFunction(3.4093+0.1066*x+0.0166*y-0.0427*z),HL,HB,HW);

SL=[LMt LPt LDt LPi LMi LDi LPm LMm LDm LPr LMr LDr LPl LMl LDl];

if G=='l'
    FF="C:\Users\panos\OneDrive - University of Warwick\PhD\Thesis\Segment length validation\";
else
    FF="C:\Users\u1857308\OneDrive - University of Warwick\PhD\Thesis\Segment length validation\";
end

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
mm=stem(L(:,j),BB(:,j),'-o');
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
    saveas(h,FF+"\Stem\segment "+f+" "+g,'jpg');
    saveas(h,FF+"\Stem\segment "+f+" "+g,'fig');
    
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
     saveas(h,FF+"\Differences\difference "+f+" "+g,'jpg');
    saveas(h,FF+"\Differences\difference "+f+" "+g,'fig');
    
    
    [H(j),p(j)]=swtest(DB(:,j)-DSL(:,j),a);
    if (H(j)==0 && p(j)>a)
        m=(DB(:,j)+DSL(:,j))/2;
        d=DB(:,j)-DSL(:,j);
        [~,pt]=ttest(DB(:,j),DSL(:,j));
            SD_B(j)=100*1.96*std(d)*2/(mean(DB(:,j))+mean(DSL(:,j)));
            s=std(d);
        t=2.262; %Value taken from the two tailed Students t-distribution table (0.05 column) from probability and statistics for engineers with 9 dof
           n=length(DB(:,j));
           b=mean(d);
           %Confidence intervals of the bias
           dci_u(j)=b+t*sqrt(s^2/n);
           dci_l(j)=b-t*sqrt(s^2/n);
           
           %Estimate the LOA
           loa_u(j)=b+1.96*s;
           loa_l(j)=b-1.96*s;
           
           %Estimate the CI of the upper LOA
           loaci_uu=loa_u(j)+t*sqrt(3*s^2/n);
           loaci_lu=loa_u(j)-t*sqrt(3*s^2/n);
           
           %Estimate the CI for the low LOA
           loaci_ul=loa_l(j)+t*sqrt(3*s^2/n);
           loaci_ll=loa_l(j)-t*sqrt(3*s^2/n);
           
           fig=figure;
           %Scatter Plot
           scatter(m,d);
           hold on;
           %Bias Line
           plot(m,b*ones(1,length(m)),'-','LineWidth',1.5);
           hold on
           %Line of Equality
           plot(m,zeros(1,length(m)),'-','LineWidth',1.5);
           hold on
           %Upper LOA
           plot(m,loa_u(j).*ones(1,length(m)),'-.','LineWidth',1.5);
           hold on
           %Lower LOA
           plot(m,loa_l(j).*ones(1,length(m)),'-.','LineWidth',1.5);
           hold on
           %Vertical line at minimum value
           plot((min(m)).*ones(1,100),linspace(dci_l(j),dci_u(j),100),'--k','LineWidth',1.5);
           hold on
           %Horizontal line at upper C.I. of bias
           
           a1=min(m)-(max(m)-min(m))*0.1;
           b1=min(m)+(max(m)-min(m))*0.1;
          
           
           plot([a1 b1],[dci_u(j) dci_u(j)],'-k','LineWidth',1.5);
           hold on;
           %Horizontal line at lower C.I. of bias
           plot([a1 b1],[dci_l(j) dci_l(j)],'-k','LineWidth',1.5);
           hold off;
           xlabel("(\Delta L_{Buchholz,Data} +\Delta L_{This study,data})/2 (mm)");
           ylabel("\Delta L_{Buchholz,Data} -\Delta L_{This study,data}(mm)");
          
            axis([a1-0.05*a1 max(m)+0.1*max(m) loa_l(j)-0.1*loa_l(j) loa_u(j)+0.1*loa_u(j)]);
           if ((loa_l(j)-0.1*loa_l(j)>0.5)|| (loa_u(j)+0.1*loa_u(j)<-0.5))
               legend("Data","Bias","","LOA upper","LOA lower","Bias C.I.",'Location','SouthEast','FontSize',8);
           else
               legend("Data","Bias","Line of Equality","LOA upper","LOA lower","Bias C.I.",'Location','SouthEast','FontSize',8);
           end
           s={"Bias: "+b+" mm","Paired t-test p-value: "+pt,"Shapiro-Wilk test p-value: "+p(j)};
           annotation('textbox', [0.3, 0.8, 0.1, 0.1],'String',s,'FitBoxToText','on')
           
     title({"Bland-Altman plot of differences between participant data and scaling functions ",...
         "from this study and Buchholz for "+f+" finger "+g+ " segment"});
     ax=gca;
     ax.FontSize=8;   
     saveas(fig,FF+"\Bland-Altman\Bland-Altman "+f+" "+g,'jpg');
        saveas(fig,FF+"\Bland-Altman\Bland-Altman "+f+" "+g,'fig');
    
    end

end


%Segment breadth validation

%Buccholz equations for segent breadth for number 2,6,10 for index to
%little and 2,6,9 for thumb. These are taken from "An ellipsoidal representaion of human hand anthropometry"
%All functions are divided by 2 since my scaling functions are for
%measuring radius and not diameter.

syms x

RMtB=feval(matlabFunction(0.225*x/2),HB);
RPtB=feval(matlabFunction(0.216*x/2),HB);
RDtB=feval(matlabFunction(0.282*x/2),HB);
RPiB=feval(matlabFunction(0.185*x/2),HB);
RMiB=feval(matlabFunction(0.206*x/2),HB);
RDiB=feval(matlabFunction(0.209*x/2),HB);
RPmB=feval(matlabFunction(0.189*x/2),HB);
RMmB=feval(matlabFunction(0.2*x/2),HB);
RDmB=feval(matlabFunction(0.202*x/2),HB);
RPrB=feval(matlabFunction(0.181*x/2),HB);
RMrB=feval(matlabFunction(0.186*x/2),HB);
RDrB=feval(matlabFunction(0.194*x/2),HB);
RPlB=feval(matlabFunction(0.167*x/2),HB);
RMlB=feval(matlabFunction(0.174*x/2),HB);
RDlB=feval(matlabFunction(0.183*x/2),HB);

BBr=[RMtB RPtB RDtB RPiB RMiB RDiB RPmB RMmB RDmB RPrB RMrB RDrB RPlB RMlB RDlB];

%My scaling equations for segment radius.

syms x y z

RMt=feval(matlabFunction(13.251+0.15*x),HW);
RPt=feval(matlabFunction(1.327+0.022*x+0.0331*y+0.114*z),PL,HB,HW);
RDt=feval(matlabFunction(-1.376+0.0148*x+0.075*y+0.066*z),HL,HB,HW);
RPi=feval(matlabFunction(5.668+0.126*x),HW);
RMi=feval(matlabFunction(0.552+0.024*x+0.045*y+0.006*z),HL,HB,HW);
RDi=feval(matlabFunction(-2.172+0.035*x+0.05*y-0.013*z),HL,HB,HW);
RPm=feval(matlabFunction(3.163+0.015*x-0.009*y+0.123*z),HL,HB,HW);
RMm=feval(matlabFunction(-0.414+0.038*x+0.006*y+0.048*z),HL,HB,HW);
RDm=feval(matlabFunction(5.172+0.086*x),HW);
RPr=feval(matlabFunction(3.419+0.016*x-0.015*y+0.165*z),PL,HB,HW);
RMr=feval(matlabFunction(0.368+0.045*x+0.029*y+0.035*z),PL,HB,HW);
RDr=feval(matlabFunction(-1.524+0.0259*x+0.043*y+0.016*z),HL,HB,HW);
RPl=feval(matlabFunction(-3.576+0.041*x+0.057*y+0.107*z),PL,HB,HW);
RMl=feval(matlabFunction(1.346+0.011*x+0.022*y+0.064*z),HL,HB,HW);
RDl=feval(matlabFunction(0.118+0.004*x+0.068*y+0.001*z),HL,HB,HW);

SR=[RMt RPt RDt RPi RMi RDi RPm RMm RDm RPr RMr RDr RPl RMl RDl];


if G=='l'
    FF1="C:\Users\panos\OneDrive - University of Warwick\PhD\Thesis\Segment radius validation\";
else
    FF1="C:\Users\u1857308\OneDrive - University of Warwick\PhD\Thesis\Segment radius validation\";
end

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
mm=stem(R(:,j),BBr(:,j),'-o');
mm.LineWidth=2;
        grid
        hold on;
        mm1=stem(R(:,j),SR(:,j),'o',':');
        mm1.LineWidth=2;
        hold on;
       grid
        xlabel("Participant segment radius data (mm)");
        ylabel("Predicted segment radius (mm)");
        title({'Graph between the predicted segment radius values of '+f+' finger',...
                ' '+g+' segment from this study, the Buchholz scaling function and the participant data'});
        legend("Buchholz function","Scaling function from this study",'Location','southeast');
        ax=gca;
        ax.FontSize=8;
    hold off;
    saveas(h,FF1+"\Stem\segment "+f+" "+g,'jpg');
    saveas(h,FF1+"\Stem\segment "+f+" "+g,'fig');
    
     clear h;
    h=figure();
    DBr(:,j)=R(:,j)-BBr(:,j);
    DSR(:,j)=R(:,j)-SR(:,j);
    plot(R(:,j),DBr(:,j),'o',R(:,j),DSR(:,j),'x');
    xlabel("Participant segment radius data (mm)");
    ylabel("Difference between participant data and predicted value (mm)");
    legend("Buchholz equation","Scaling equation from this study",'Location','southeast');
    title({'Difference between participant data and scaling equation',... 
            'for '+f+' finger '+g+' segment'});
     ax=gca;
     ax.FontSize=8;
     saveas(h,FF1+"\Differences\difference "+f+" "+g,'jpg');
    saveas(h,FF1+"\Differences\difference "+f+" "+g,'fig');
    
    
    [H(j),p(j)]=swtest(DBr(:,j)-DSR(:,j),a);
    if (H(j)==0 && p(j)>a)
        m=(DBr(:,j)+DSR(:,j))/2;
        d=DBr(:,j)-DSR(:,j);
        [~,pt]=ttest(DBr(:,j),DSR(:,j));
            SD_Br(j)=100*1.96*std(d)*2/(mean(DBr(:,j))+mean(DSR(:,j)));
            s=std(d);
        t=2.262; %Value taken from t student table 0.025 from probability and statistics for engineers with 9 dof
           n=length(DBr(:,j));
           b=mean(d);
           %Confidence intervals of the bias
           dci_u(j)=b+t*sqrt(s^2/n);
           dci_l(j)=b-t*sqrt(s^2/n);
           
           %Estimate the LOA
           loa_u(j)=b+1.96*s;
           loa_l(j)=b-1.96*s;
           
           %Estimate the CI of the upper LOA
           loaci_uu=loa_u(j)+t*sqrt(3*s^2/n);
           loaci_lu=loa_u(j)-t*sqrt(3*s^2/n);
           
           %Estimate the CI for the low LOA
           loaci_ul=loa_l(j)+t*sqrt(3*s^2/n);
           loaci_ll=loa_l(j)-t*sqrt(3*s^2/n);
           
           fig=figure;
           %Scatter Plot
           scatter(m,d);
           hold on;
           %Bias Line
           plot(m,b*ones(1,length(m)),'-','LineWidth',1.5);
           hold on
           %Line of Equality
           plot(m,zeros(1,length(m)),'-','LineWidth',1.5);
           hold on
           %Upper LOA
           plot(m,loa_u(j).*ones(1,length(m)),'-.','LineWidth',1.5);
           hold on
           %Lower LOA
           plot(m,loa_l(j).*ones(1,length(m)),'-.','LineWidth',1.5);
           hold on
           %Vertical line at minimum value
           plot((min(m)).*ones(1,100),linspace(dci_l(j),dci_u(j),100),'--k','LineWidth',1.5);
           hold on
           %Horizontal line at upper C.I. of bias
           
           a1=min(m)-(max(m)-min(m))*0.1;
           b1=min(m)+(max(m)-min(m))*0.1;
          
           
           plot([a1 b1],[dci_u(j) dci_u(j)],'-k','LineWidth',1.5);
           hold on;
           %Horizontal line at lower C.I. of bias
           plot([a1 b1],[dci_l(j) dci_l(j)],'-k','LineWidth',1.5);
           hold off;
           xlabel("(\Delta R_{Buchholz,Data} +\Delta R_{This study,data})/2 (mm)");
           ylabel("\Delta R_{Buchholz,Data} -\Delta R_{This study,data}(mm)");
          
            axis([a1-0.05*a1 max(m)+0.1*max(m) loa_l(j)-0.1*loa_l(j) loa_u(j)+0.1*loa_u(j)]);
           if ((loa_l(j)-0.1*loa_l(j)>0.5)|| (loa_u(j)+0.1*loa_u(j)<-0.5))
               legend("Data","Bias","","LOA upper","LOA lower","Bias C.I.",'Location','SouthEast','FontSize',8);
           else
               legend("Data","Bias","Line of Equality","LOA upper","LOA lower","Bias C.I.",'Location','SouthEast','FontSize',8);
           end
           s={"Bias: "+b+" mm","Paired t-test p-value: "+pt,"Shapiro-Wilk test p-value: "+p(j)};
           annotation('textbox', [0.3, 0.8, 0.1, 0.1],'String',s,'FitBoxToText','on')
           
     title({"Bland-Altman plot of differences between participant data and scaling functions ",...
         "from this study and Buchholz for "+f+" finger "+g+ " segment"});
     ax=gca;
     ax.FontSize=8;   
     saveas(fig,FF1+"\Bland-Altman\Bland-Altman "+f+" "+g,'jpg');
        saveas(fig,FF1+"\Bland-Altman\Bland-Altman "+f+" "+g,'fig');
    
    end

end
