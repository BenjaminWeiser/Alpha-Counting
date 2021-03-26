clear;
%Made by Benjamin Weiser, April 2020.

% Uses functions findranges.m, poly4fit2fun.m, vectorizefunction.m

% See 2020 Summer work placement report for more info.

%This code is looking at determing coprecipitation yields from alpha
%counting measurments. 

% Since the trace amounts of 233U and 229Th added are too dilute to
% precipitate by conventional means, it is typically coprecipitated with a
% carrier (Cu(OH)2), a substance that has a similar crystalline
% structurethat can incorporate the desired element. The coprecipitation
% yield is the percent incorporation of the actinide into the carrier

%By alpha counting we will be able to measure the concentration and thus
%the coprecipitationyields of both U and Th with just one sample.

%Since alpha particles are charged particles, they lose energy in the
%medium, therefore most of them areabsorbed by the material itself and not
%emitted out of the material. Energy deposition of the emitted
%alphaparticles for each atom in the decay chain was simulated in SRIM
%using the most prevalent alpha decayenergy.



%Set parameters
%Energy resolution of alpha counter (minimum observable energy)
Eremaining=1.5E6;

%initial energies alpha decay energies

%U232 chain
EiU233=4800E3;
EiU232=5300E3;
EiTh228=5400E3;
EiRa224=5685E3;
EiRn220=6288E3;
EiPo216=6778E3;
EiBi212=6050E3;
EiPo212=8785E3;

%Th229 chain
EiTh229=4845E3;
EiAc225=5830E3;
EiFr221=6341E3;
EiAt217=7066E3;
EiPo213=8376E3;


%Halflives (years)
TU233=1.59E5;
TU232=6.89E1;
TTh228=1.91E0;
TTh229=7932;
T =[TU233 TU232 TTh228 TTh229];


%Energy deposition of the emitted alphaparticles for each atom in the decay
%chain was simulated in SRIM using the most prevalent alpha decayenergy.
%Read files
U233 = xlsread('Activity_of_1ppt_U233.xlsx', 'A27:E126');
U232 = xlsread('activity of U232.xlsx', 'A27:E126');
Th228 = xlsread('TH228alphapenetration.xlsx', 'A27:E126');
Ra224 = xlsread('alphaRa224', 'A27:E126');
Rn220 = xlsread('alpha220Rn', 'A27:E126');
Po216 = xlsread('alpha216Po', 'A27:E126');
Bi212 = xlsread('alphaBi212', 'A27:E126');
Po212 = xlsread('alphaPo212', 'A27:E126');
%Read file for 229Th chain
Th229 = xlsread('alpha229Th', 'A27:E126');
Ac225 = xlsread('alpha225Ac', 'A27:E126');
Fr221 = xlsread('alpha221Fr', 'A27:E126');
At217 = xlsread('alpha217At', 'A27:E126');
Po213 = xlsread('alpha213Po', 'A27:E126');


%The integral of the energy deposition over the depth provides the energy
%lost by the alpha particle.Subtracting the energy lost from the initial
%energy of the alpha’s we obtain the energy distribution whichwe will
%observe in the alpha counter. Since this counter has a minimum observable
%energy of 1.5MeV wetake the distance for the alpha particles to reach that
%energy as the max range in which the particle cantraverse through the
%copper and still be counted by the detector. This allows us to determine
%the volume ofthe sample which contains the radionuclides we observe by the
%detector.

%%set arrays and find ranges (RX###=max range) 
U233= U233(1:end,1:2);%column 1 is depth(angtrom), column 2 is Eloss/depth i eV/Ang
RU233 = findranges(U233,EiU233);

U232= U232(1:end,1:2);%column 1 is depth(angtrom), column 2 is Eloss/depth i eV/Ang
RU232 = findranges(U232,EiU232);

Th228 = Th228(1:end,1:2);%column 1 is depth(angtrom), column 2 is Eloss/depth i eV/Ang
RTH228 = findranges(Th228,EiTh228);

Ra224 = Ra224(1:end,1:2);%column 1 is depth(angtrom), column 2 is Eloss/depth i eV/Ang
RRa224 = findranges(Ra224,EiRa224);

Rn220 = Rn220(1:end,1:2);%column 1 is depth(angtrom), column 2 is Eloss/depth i eV/Ang
RRn220 = findranges(Rn220,EiRn220);

Po216 = Po216(1:end,1:2);%column 1 is depth(angtrom), column 2 is Eloss/depth i eV/Ang
RPo216 = findranges(Po216,EiPo216);

Bi212 = Bi212(1:end,1:2);%column 1 is depth(angtrom), column 2 is Eloss/depth i eV/Ang
RBi212 = findranges(Bi212,EiBi212);

Po212 = Po212(1:end,1:2);%column 1 is depth(angtrom), column 2 is Eloss/depth i eV/Ang
RPo212 = findranges(Po212,EiPo212);
%%%%% Th229 chain
Th229= Th229(1:end,1:2);%column 1 is depth(angtrom), column 2 is Eloss/depth i eV/Ang
RTh229 = findranges(Th229,EiTh229);

Ac225= Ac225(1:end,1:2);%column 1 is depth(angtrom), column 2 is Eloss/depth i eV/Ang
RAc225 = findranges(Ac225,EiAc225);

Fr221= Fr221(1:end,1:2);%column 1 is depth(angtrom), column 2 is Eloss/depth i eV/Ang
RFr221 = findranges(Fr221,EiFr221);

At217= At217(1:end,1:2);%column 1 is depth(angtrom), column 2 is Eloss/depth i eV/Ang
RAt217 = findranges(At217,EiAt217);

Po213= Po213(1:end,1:2);%column 1 is depth(angtrom), column 2 is Eloss/depth i eV/Ang
RPo213 = findranges(Po213,EiPo213);


%puts ranges in array
RangesU232=[RU232(1,3) RTH228(1,3) RRa224(1,3) RRn220(1,3) RPo216(1,3) RBi212(1,3) RPo212(1,3)];
RangesTh229=[RTh229(1,3) RAc225(1,3) RFr221(1,3) RAt217(1,3) RPo213(1,3)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This section is to calculate observed activity. The distance between the
%sample and the detector can be approximated to be 0 since the distance is
%soshort meaning that the activity observed by the detector is one half of
%the total activity.Using the surface area of the target multiplied by the
%max range of the alpha particles allows us to get theactivity observed for
%each radionuclide in the chain. constants

% Theseactivities are then divided by their respective maximum alpha range
% to give activity per target depth dA/dx which is constant. The depth over
% energy data was then fitted with a 4th order polynomial. The
% derivative was then calculated yieldingdxdE. The result of multiplying the
% two is a function dA/dE(E)for each actinide.Summing each of the function
% gives the prediction for what we will observe by the detector.
% Multiplying each of the daughters of232U by a yield factor allows us to
% see the difference that the 228Th coprecipitationyield will give.


M_Cu=63.55; %g/mol
rho_Cu=8.92; %g/cm^3
avagadros=6.02214076E23; %atom/mol

%Concentrations placed in samples
concentrationU232=0.1E-12; %1ppt, concentration of the actinide g/g
concentrationTh229=0.1E-12;
M_U232=232; %mass of that actinide
M_Th229=229;


%Defining geometry of target and detector
dt=1; %diameter of target (cm)
d=0.1; %distance between target and detector in (cm)
r=15; %radius of the detector in (cm)




%Finding activityy using initial activity but differentranges for each
%actinide due to difference in initial alpha energies

%index of i choses actinide with starting concentration. i picks index in
%ranges and in halflife
%calculating actvity per volume using T(i) and concentration where i picks
%the actinide. i=2 is for U232, i=4 is for Th229. 

Cu_per_V=8.92/63.55*avagadros; %#Cu/cm^3
Actinide_per_V_U232=Cu_per_V*concentrationU232*M_Cu/M_U232; %#U/cm^3
Actinide_per_V_Th229=Cu_per_V*concentrationTh229*M_Cu/M_Th229;

Activity_per_V_U232=Actinide_per_V_U232*log(2)/T(2)/8760; %in alpha/h/cm^3
Activity_per_V_bq=Activity_per_V_U232/60/60;% in bq/cm^3
Activity_per_V_Th229=Actinide_per_V_Th229*log(2)/T(4)/8760; %in alpha/h/cm^3

Activity_from_target_U232=Activity_per_V_U232*(dt/2)^2*pi.*RangesU232/1E8; %uses actvity of mother and ranges of each daughter
Activity_from_target_bq=Activity_per_V_bq*(dt/2)^2*pi.*RangesU232;
Activity_from_target_Th229=Activity_per_V_Th229*(dt/2)^2*pi.*RangesTh229/1E8;

%U232
Activity_detector_U232=1/2*(-d/sqrt(d^2+r^2)+1)*Activity_from_target_U232; %alpha/h
Activity_detector_U232(1,6)=36/100*Activity_detector_U232(1,6); %Bi alpha only 36% of the time
Activity_detector_U232(1,7)=64/100*Activity_detector_U232(1,7); %Po from beta of Bi 64% of the time
%Th229
Activity_detector_Th229=1/2*(-d/sqrt(d^2+r^2)+1)*Activity_from_target_Th229; 

%fit x versus E graph
 
%ranges are in A. Convert A to um
Constant_U232=Activity_detector_U232./(RangesU232*1E-4); %alpha/hour/um
Constant_Th229=Activity_detector_Th229./(RangesTh229*1E-4);

%poly4fit2fun converts units to MeV and to micrometer. It fits the X versus E
%graph and returns the coefficients
Fun_U232 = poly4fit2fun(RU232).*Constant_U232(1,1);  %coefficients for graph of x versus E
Fun_TH228= poly4fit2fun(RTH228).*Constant_U232(1,2);
Fun_Ra224= poly4fit2fun(RRa224).*Constant_U232(1,3);
Fun_Rn220= poly4fit2fun(RRn220).*Constant_U232(1,4);
Fun_Po216= poly4fit2fun(RPo216).*Constant_U232(1,5);
Fun_Bi212= poly4fit2fun(RBi212).*Constant_U232(1,6);
Fun_Po212= poly4fit2fun(RPo212).*Constant_U232(1,7);

%Th229 chain
Fun_Th229= poly4fit2fun(RPo212).*Constant_Th229(1,1);
Fun_Ac225= poly4fit2fun(RPo212).*Constant_Th229(1,2);
Fun_Fr221= poly4fit2fun(RPo212).*Constant_Th229(1,3);
Fun_At217= poly4fit2fun(RPo212).*Constant_Th229(1,4);
Fun_Po213= poly4fit2fun(RPo212).*Constant_Th229(1,5);

%take derivitive of x/E to get dX/dE
%funtion give units of alpha/hour as a funciton of E

%index of i choses actinide with starting concentration. i picks index in
%ranges and in halflife
%first activity calculations used
%i=1; 
%{
for E=1.5:0.1:10

    A=0;
    AU232 = -(4*Fun_U232(1,1)*E.^3 + 3*Fun_U232(1,2)*E.^2 + 2*Fun_U232(1,3)*E + Fun_U232(1,4));
    ATH228 = -(4*Fun_TH228(1,1)*E.^3 + 3*Fun_TH228(1,2)*E.^2 + 2*Fun_TH228(1,3)*E + Fun_TH228(1,4));
    ARa224 = -(4*Fun_Ra224(1,1)*E.^3 + 3*Fun_Ra224(1,2)*E.^2 + 2*Fun_Ra224(1,3)*E + Fun_Ra224(1,4));
    ARn220 = -(4*Fun_Rn220(1,1)*E.^3 + 3*Fun_Rn220(1,2)*E.^2 + 2*Fun_Rn220(1,3)*E + Fun_Rn220(1,4));
    ABi212 = -(4*Fun_Bi212(1,1)*E.^3 + 3*Fun_Bi212(1,2)*E.^2 + 2*Fun_Bi212(1,3)*E + Fun_Bi212(1,4));
    APo212 = -(4*Fun_Po212(1,1)*E.^3 + 3*Fun_Po212(1,2)*E.^2 + 2*Fun_Po212(1,3)*E + Fun_Po212(1,4));
    if E<5.3
    A=A+AU232;
    end
    if E<5.4
    A=A+ATH228;
    end
    if E<5.685
    A=A+ARa224;
    end
    if E<6.288
    A=A+ARn220;
    end
    if E<6.050
    A=A+ABi212;
    end
    if E<8.785
    A=A+APo212;
    end
 Atotal(i,1)= E;
 Atotal(i,2)= A;
 i=i+1;
end
%}
 %derivative give alpha per hour per MeV
%For A=dN/dE. where Yield is %yield. doing gor 20% 60% and 100% of Th to U.
%Also makes specific array for only values greater that E initial of U232
n=1;

for Yield=0.20:0.40:1
i=1; 
m=1;                                                                                                                                                                                                            
for E=1.5:0.0001:10
    A=0;
    B=0;
    AU232 = -(4*Fun_U232(1,1)*E.^3 + 3*Fun_U232(1,2)*E.^2 + 2*Fun_U232(1,3)*E + Fun_U232(1,4));
    ATH228 = -Yield*(4*Fun_TH228(1,1)*E.^3 + 3*Fun_TH228(1,2)*E.^2 + 2*Fun_TH228(1,3)*E + Fun_TH228(1,4));
    ARa224 = -Yield*(4*Fun_Ra224(1,1)*E.^3 + 3*Fun_Ra224(1,2)*E.^2 + 2*Fun_Ra224(1,3)*E + Fun_Ra224(1,4));
    ARn220 = -Yield*(4*Fun_Rn220(1,1)*E.^3 + 3*Fun_Rn220(1,2)*E.^2 + 2*Fun_Rn220(1,3)*E + Fun_Rn220(1,4));
    APo216 = -Yield*(4*Fun_Po216(1,1)*E.^3 + 3*Fun_Po216(1,2)*E.^2 + 2*Fun_Po216(1,3)*E + Fun_Po216(1,4));
    ABi212 = -Yield*(4*Fun_Bi212(1,1)*E.^3 + 3*Fun_Bi212(1,2)*E.^2 + 2*Fun_Bi212(1,3)*E + Fun_Bi212(1,4));
    APo212 = -Yield*(4*Fun_Po212(1,1)*E.^3 + 3*Fun_Po212(1,2)*E.^2 + 2*Fun_Po212(1,3)*E + Fun_Po212(1,4));
   %Th229 chain
    ATh229 = -Yield*(4*Fun_Th229(1,1)*E.^3 + 3*Fun_Th229(1,2)*E.^2 + 2*Fun_Th229(1,3)*E + Fun_Th229(1,4));
    AAc225 = -Yield*(4*Fun_Ac225(1,1)*E.^3 + 3*Fun_Ac225(1,2)*E.^2 + 2*Fun_Ac225(1,3)*E + Fun_Ac225(1,4));
    AFr221 = -Yield*(4*Fun_Fr221(1,1)*E.^3 + 3*Fun_Fr221(1,2)*E.^2 + 2*Fun_Fr221(1,3)*E + Fun_Fr221(1,4));
    AAt217 = -Yield*(4*Fun_At217(1,1)*E.^3 + 3*Fun_At217(1,2)*E.^2 + 2*Fun_At217(1,3)*E + Fun_At217(1,4));
    APo213 = -Yield*(4*Fun_Po213(1,1)*E.^3 + 3*Fun_Po213(1,2)*E.^2 + 2*Fun_Po213(1,3)*E + Fun_Po213(1,4));
    if E<5.3
    A=A+AU232;
    B=B+AU232;
    end
    if E<5.4
    A=A+ATH228;
    B=B+ATH228;
    end
    if E<5.685
    A=A+ARa224;
    B=B+ARa224;
    end
    if E<6.288
    A=A+ARn220;
    B=B+ARn220;
    end
    if E<6.778
    A=A+APo216;
    B=B+APo216;
    end
    if E<6.050
    A=A+ABi212;
    B=B+ABi212;
    end
    if E<8.785
    A=A+APo212;
    B=B+APo212;
    end
    %Th229 chain
    if E<4.845
    A=A+ATh229;
    end
    if E<5.830
    A=A+AAc225;
    end
    if E<6.341
    A=A+AFr221;
    end
    if E<7.0669
    A=A+AAt217;
    end
    if E<8.376
    A=A+APo213;
    end
 AtotalYield(i,n)= E;
 AtotalYield(i,n+1)= A;
 AU232Chain(i,n)= E;
 AU232Chain(i,n+1)= B;
 if E>7.5
     Agreaterthan75(m,n)= E;
     Agreaterthan75(m,n+1)= A;
     m=m+1;
 end
 i=i+1;
end
n=n+2;
end
  
   FU232 = @(E) -(4*Fun_U232(1,1)*E.^3 + 3*Fun_U232(1,2)*E.^2 + 2*Fun_U232(1,3)*E + Fun_U232(1,4));
   FTh228 = @(E)-(4*Fun_TH228(1,1)*E.^3 + 3*Fun_TH228(1,2)*E.^2 + 2*Fun_TH228(1,3)*E + Fun_TH228(1,4));
   % FRa224 = @(E)-(4*Fun_Ra224(1,1)*E.^3 + 3*Fun_Ra224(1,2)*E.^2 + 2*Fun_Ra224(1,3)*E + Fun_Ra224(1,4));
   % FRn220 = @(E)-(4*Fun_Rn220(1,1)*E.^3 + 3*Fun_Rn220(1,2)*E.^2 + 2*Fun_Rn220(1,3)*E + Fun_Rn220(1,4));
   % FPo216 = @(E)-(4*Fun_Po216(1,1)*E.^3 + 3*Fun_Po216(1,2)*E.^2 + 2*Fun_Po216(1,3)*E + Fun_Po216(1,4));
   % FBi212 = @(E)-(4*Fun_Bi212(1,1)*E.^3 + 3*Fun_Bi212(1,2)*E.^2 + 2*Fun_Bi212(1,3)*E + Fun_Bi212(1,4));
   % FPo212 = @(E)-(4*Fun_Po212(1,1)*E.^3 + 3*Fun_Po212(1,2)*E.^2 + 2*Fun_Po212(1,3)*E + Fun_Po212(1,4));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% data was used to calculate two constants to determine the specific alpha
% particles coming from 232U and 228Th from the entire decay chain alpha
% particles. The alpha energy of the 232U nuclide is 5.3MeV andis the lowest
% of the energies in the decay chain. From the alpha counting measurement,
% we will know thetotal rate of alpha particles reaching the detector,?
%Total, and the alpha per hour which hit the detector withgreater energy
% than 7.5MeV,? E>7.5MeV.
   
%vectorizes U232 function     
VU232=vectorizefunction(FU232, 1.5, 0.0001, 5.3);  
VU232(38001,1)=5.3;
VU232(38001,2)=0;
TotalU232=trapz(VU232(:,1),VU232(:,2)); 

VTh228=vectorizefunction(FTh228, 1.5, 0.0001, 5.4);  
VTh228(38001,1)=5.4;
VTh228(38001,2)=0;
TotalTh228=trapz(VTh228(:,1),VTh228(:,2)); 

Sumofactivities=sum(Activity_detector_U232)+ sum(Activity_detector_Th229) %alpha/h, gives 69.1724 for 1ppt U232 chain, 
%total with 0.1 ppt U232 and 1ppt Th229 is 8.3989. U232 contribution is
%6.917. Sum activities for 1ppt Th229 chain is 1.4816
n=1;
i=1;
while i<=6
Totalalphas(1,n)=trapz(AtotalYield(:,i),AtotalYield(:,i+1)); %gives 68.4356
TotalgreaterthanU232(1,n)=trapz(Agreaterthan75(:,i),Agreaterthan75(:,i+1));
n=n+1;
i=i+2;
end
%alphasgreaterthan5.3eV=trapz(AtotalYield(:,5),AtotalYield(:,6))

%constant=th228chain/greaterthan5.3
U232Constant=(Totalalphas-TotalU232)./TotalgreaterthanU232; %alpha/alpha
Th228Constantyield20=0.20*TotalTh228./TotalgreaterthanU232(1,1);
Th228Constantyield60=0.60*TotalTh228./TotalgreaterthanU232(1,2);
Th228Constantyield100=TotalTh228./TotalgreaterthanU232(1,3);


%Check to confirm code!!!
%calculating concentration opposite wayy. Should get 1ppt
alphaU232=Totalalphas-U232Constant.*TotalgreaterthanU232
C=alphaU232*M_U232*TU232*(8760)/(log(2)*rho_Cu*avagadros*(dt/2)^2*pi*(RangesU232(1,1)/1E8))




figure
grid on; hold on;
%plot(RU232(:,2),RU232(:,1),'x', RTH228(:,2), RTH228(:,1),'*', RRa224(:,2), RRa224(:,1),'*', RRn220(:,2), RRn220(:,1),'*', RBi212(:,2), RBi212(:,1),'*', RPo212(:,2), RPo212(:,1),'*')
plot(AtotalYield(:,1),AtotalYield(:,2),AtotalYield(:,3),AtotalYield(:,4),AtotalYield(:,5),AtotalYield(:,6))
 %fplot(FU232)
 %fplot(FTH228)
 %fplot(FRa224) 
 %fplot(FRn220) 
 %fplot(FBi212) 
 %fplot(FPo212)
 xlim([1.5 10])
 %ylim([0 60])
 ylabel('Energy distribution of the alpha particles ($\frac{1}{h MeV}$)', 'interpreter', 'latex', 'fontsize', 10);
xlabel('Energy of detected alphas ($MeV$)', 'interpreter', 'latex','fontsize', 10);
title({'Change in activity versus Energy of Emitted Alphas','for copper sample with 0.1ppt 232U, 0.1ppt Th229','with different coprecipiation yields of Th'}, 'interpreter', 'latex', 'fontsize', 10);

%ledge = legend("Total Activity", "U232","Th228","Ra224","Rn220","Bi212","Po212");
ledge = legend("20% Th Yield","60% Th Yield","100% Th Yield");
ledge.Location = 'northeast';

figure
grid on; hold on;
plot(AtotalYield(:,5),AtotalYield(:,6),AU232Chain(:,5),AU232Chain(:,6))
 %fplot(FU232)
 %fplot(FTH228)
 %fplot(FRa224) 
 %fplot(FRn220) 
 %fplot(FBi212) 
 %fplot(FPo212)
 xlim([1.5 10])
 %ylim([0 60])
 ylabel('Energy distribution of the alpha particles ($\frac{1}{h MeV}$)', 'interpreter', 'latex', 'fontsize', 10);
xlabel('Energy of detected alphas ($MeV$)', 'interpreter', 'latex','fontsize', 10);
title({'Change in activity versus energy of emitted alphas','for copper sample with 0.1ppt 232U','with and without addition of 0.1ppt Th229'}, 'interpreter', 'latex', 'fontsize', 10);

%ledge = legend("Total Activity", "U232","Th228","Ra224","Rn220","Bi212","Po212");
ledge = legend("With addtion of 0.1ppt Th229","Only 232U Chain");
ledge.Location = 'northeast';

