% function [dx_dt]= eqn1(t,x)
% %a function which returns a rate of change vector
% % global AlpaIme1 AlpaIme2 AlpaCln3 AlpaRpd3Sin3 BetaIme1 BetaIme2 BetaCln3 n1 n2 n3 n4 n5 n6 n7 n8 KIme11 KIme22 KN2 KIme12 KRpd3 GammaPIme2 KPIme2 dIme2 dIme1 KN2Cln3 GammaPCdk1 KPIme2Cdk1 dCln3 dRpd3  N2 KCln3Ime1 n9 GammaPIme1Cln3 n10 KPIme1Cln3 KRpd3Ime1 n12 n18 xo yo po zo no ;
%   parameters
% % % Nitrogen depletion as a constant high value
% % dx_dt(1) = (AlpaIme1/xo)+(BetaIme1/xo)*(((x(1))^n1)/((KIme11/xo)^n1+(x(1))^n1))*(KN2^(n2)/(KN2^n2+(N2*no)^n2))*((KCln3Ime1^(n9))/(KCln3Ime1^n9+(x(3)*zo)^n9))-(GammaPIme1Cln3/xo)*(((x(3))^n10)/((KPIme1Cln3/zo)^n10+(x(3))^n10))-(dIme1*x(2)*yo*x(1))^(1) ;
% % dx_dt(2) = (AlpaIme2/yo)+(BetaIme2/yo)*((((x(2))^n3)/((KIme22/yo)^n3+(x(2))^n3))*(((x(1))^n4)/((KIme12/xo)^n4+(x(1))^n4))*(((KRpd3)^n5)/((KRpd3)^n5+(x(4)*po)^n5)))-(GammaPIme2/yo)*(((x(3))^n6)/((KPIme2/zo)^n6+(x(3))^n6))- dIme2*x(2);
% % %dx_dt(2) = (AlpaIme2/yo)+(BetaIme2/yo)*((((x(1))^n4)/((KIme12/xo)^n4+(x(1))^n4))*(((KRpd3)^n5)/((KRpd3)^n5+(x(4)*po)^n5)))-(GammaPIme2/yo)*(((x(3))^n6)/((KPIme2/zo)^n6+(x(3))^n6))- dIme2*x(2);
% % dx_dt(3) = (AlpaCln3/zo)+(BetaCln3/zo)*((N2))^(n7)/((KN2Cln3/no)^n7+(N2)^n7)- (GammaPCdk1/zo)*(((x(2))^n8)/((KPIme2Cdk1/yo)^n8+(x(2))^n8))-dCln3*x(3);
% % dx_dt(4) = (AlpaRpd3Sin3/po)-dRpd3*((x(1)*xo)/((N2*no*po)));
% 
% % Nitrogen depletion as a constant high value
% % dx_dt(1) = (1/1)*(AlpaIme1/xo)+(1/1)*(BetaIme1/xo)*(((x(1))^n1)/((KIme11/xo)^n1+(x(1))^n1))*(KN2^(n2)/(KN2^n2+(N2*no)^n2))*((KCln3Ime1^(n9))/(KCln3Ime1^n9+(x(3)*zo)^n9))*((KRpd3Ime1^(n12))/(KRpd3Ime1^n12+((n18-(n18/(1+exp(t-6))))*po)^n12))-(dIme1*x(1))^(1) ;
% dx_dt(1) = (3/3)*(AlpaIme1/xo)+(3/3)*(BetaIme1/xo)*(((x(1))^n1)/((KIme11/xo)^n1+(x(1))^n1))*(KN2^(n2)/(KN2^n2+(N2*no)^n2))*((KCln3Ime1^(n9))/(KCln3Ime1^n9+(x(3)*zo)^n9))*((KRpd3Ime1^(n12))/(KRpd3Ime1^n12+((n18-(n18/(1+exp(t-6))))*po)^n12))-(GammaPIme1Cln3/xo)*(((x(3))^n10)/((KPIme1Cln3/zo)^n10+(x(3))^n10))-(dIme1*x(2)*yo*x(1))^(1) ;
% % dx_dt(1) = (3/3)*(AlpaIme1/xo)+(3/3)*(BetaIme1/xo)*(((x(1))^n1)/((KIme11/xo)^n1+(x(1))^n1))*(KN2^(n2)/(KN2^n2+(N2*no)^n2))*((KCln3Ime1^(n9))/(KCln3Ime1^n9+(x(3)*zo)^n9))*((KRpd3Ime1^(n12))/(KRpd3Ime1^n12+((n18-(n18/(1+exp(t-6))))*po)^n12))-(GammaPIme1Cln3/xo)*(((x(3))^n10)/((KPIme1Cln3/zo)^n10+(x(3))^n10))-(dIme1*x(2)*yo*x(1))^(1) ;
% 
% dx_dt(2) = (AlpaIme2/yo)+(BetaIme2/yo)*((((x(2))^n3)/((KIme22/yo)^n3+(x(2))^n3))*(((x(1))^n4)/((KIme12/xo)^n4+(x(1))^n4))*(((KRpd3)^n5)/((KRpd3)^n5+(x(4)*po)^n5)))-(GammaPIme2/yo)*(((x(3))^n6)/((KPIme2/zo)^n6+(x(3))^n6))-dIme2*x(2);
% % dx_dt(2) = (AlpaIme2/yo)+(BetaIme2/yo)*((((x(2))^n3)/((KIme22/yo)^n3+(x(2))^n3))*(((x(1))^n4)/((KIme12/xo)^n4+(x(1))^n4))*(((KRpd3)^n5)/((KRpd3)^n5+(x(4)*po)^n5)))-dIme2*x(2);
% 
% dx_dt(3) = (AlpaCln3/zo)+(BetaCln3/zo)*((N2))^(n7)/((KN2Cln3/no)^n7+(N2)^n7)-(GammaPCdk1/zo)*(((x(2))^n8)/((KPIme2Cdk1/yo)^n8+(x(2))^n8))-dCln3*x(3);
% % dx_dt(3) = (AlpaIme2/yo)+(BetaIme2/yo)*((((x(2))^n3)/((KIme22/yo)^n3+(x(2))^n3))*(((x(1))^n4)/((KIme12/xo)^n4+(x(1))^n4))*(((KRpd3)^n5)/((KRpd3)^n5+(x(4)*po)^n5)))- (GammaPCdk1/zo)*(((x(2))^n8)/((KPIme2Cdk1/yo)^n8+(x(2))^n8))-dCln3*x(3);
% % dx_dt(3) = (AlpaCln3/zo)+(BetaCln3/zo)*((N2))^(n7)/((KN2Cln3/no)^n7+(N2)^n7)-dCln3*x(3);
% 
% % dx_dt(4) = (AlpaRpd3Sin3/po)-dRpd3*((x(1)*xo)*1/po)*((410/(0.01+(exp((1.6*t))))+0.06));
% % dx_dt(4) = (AlpaRpd3Sin3/po)-dRpd3*((x(1)*xo)*1/po)*((1/((N2*no)))*(296/(0.01+(exp((1.6*t))))+0.06));
% dx_dt(4) = (AlpaRpd3Sin3/po)-dRpd3*((x(1)*xo)*1/po)*((1/((2*N2*no)))*((296)/(0.01+(exp((1.6*t))))+0.06)); % Changed to get the Nitrogen in 0 to 1 range.
% 
% % dx_dt(1) = (AlpaIme1/xo)+(BetaIme1/xo)*(((x(1))^n1)/( (KIme11/xo)^n1+(x(1))^n1))*(KN2^(n2)/(KN2^n2+((1/(8+exp(t*(0.1)-0.01)))*no)^n2))*((KCln3Ime1^(n9))/(KCln3Ime1^n9+(x(3)*zo)^n9))-(GammaPIme1Cln3/xo)*(((x(3))^n10)/((KPIme1Cln3/zo)^n10+(x(3))^n10))-(dIme1*x(2)*yo*x(1))^(1) ;
% % dx_dt(2) = (AlpaIme2/yo)+(BetaIme2/yo)*((((x(2))^n3)/((KIme22/yo)^n3+(x(2))^n3))*(((x(1))^n4)/((KIme12/xo)^n4+(x(1))^n4))*(((KRpd3)^n5)/((KRpd3)^n5+(x(4)*po)^n5)))-(GammaPIme2/yo)*(((x(3))^n6)/((KPIme2/zo)^n6+(x(3))^n6))- dIme2*x(2);
% % %dx_dt(2) = (AlpaIme2/yo)+(BetaIme2/yo)*((((x(1))^n4)/((KIme12/xo)^n4+(x(1))^n4))*(((KRpd3)^n5)/((KRpd3)^n5+(x(4)*po)^n5)))-(GammaPIme2/yo)*(((x(3))^n6)/((KPIme2/zo)^n6+(x(3))^n6))- dIme2*x(2);
% % 
% % dx_dt(3) =(AlpaCln3/zo)+(BetaCln3/zo)*(((1/(8+exp(t*(0.1)-0.01)))))^(n7)/((KN2Cln3/no)^n7+((1/(8+exp(t*(0.1)-0.01))))^n7)- (GammaPCdk1/zo)*(((x(2))^n8)/((KPIme2Cdk1/yo)^n8+(x(2))^n8))-dCln3*x(3);
% % dx_dt(4)=(AlpaRpd3Sin3/po)-dRpd3*((x(1)*xo)/(((1/(8+exp(t*(0.1)-0.01))))*no*po));
% 
% % dx_dt(1) = (AlpaIme1/xo)+(BetaIme1/xo)*(((x(1))^n1)/((KIme11/xo)^n1+(x(1))^n1))*(KN2^(n2)/(KN2^n2+((0.3/(0.05+(exp((2*t))))+0.05)*no)^n2))*((KCln3Ime1^(n9))/(KCln3Ime1^n9+(x(3)*zo)^n9))-(GammaPIme1Cln3/xo)*(((x(3))^n10)/((KPIme1Cln3/zo)^n10+(x(3))^n10))-(dIme1*x(2)*yo*x(1))^(1) ;
% % dx_dt(2) = (AlpaIme2/yo)+(BetaIme2/yo)*((((x(2))^n3)/((KIme22/yo)^n3+(x(2))^n3))*(((x(1))^n4)/((KIme12/xo)^n4+(x(1))^n4))*(((KRpd3)^n5)/((KRpd3)^n5+(x(4)*po)^n5)))-(GammaPIme2/yo)*(((x(3))^n6)/((KPIme2/zo)^n6+(x(3))^n6))- dIme2*x(2);
% % %dx_dt(2) = (AlpaIme2/yo)+(BetaIme2/yo)*((((x(1))^n4)/((KIme12/xo)^n4+(x(1))^n4))*(((KRpd3)^n5)/((KRpd3)^n5+(x(4)*po)^n5)))-(GammaPIme2/yo)*(((x(3))^n6)/((KPIme2/zo)^n6+(x(3))^n6))- dIme2*x(2);
% % 
% % dx_dt(3) = (AlpaCln3/zo)+(BetaCln3/zo)*(((0.3/(0.05+(exp((2*t))))+0.05)))^(n7)/((KN2Cln3/no)^n7+(0.3/(0.05+(exp((2*t))))+0.05)^n7)- (GammaPCdk1/zo)*(((x(2))^n8)/((KPIme2Cdk1/yo)^n8+(x(2))^n8))-dCln3*x(3);
% % dx_dt(4) = (AlpaRpd3Sin3/po)-dRpd3*(((x(1)*1))*xo)*1/(((0.3/(0.05+(exp((2*t))))+0.05))*no*po);
% 
% % Transpose dx_dt so it is a column vector
% dx_dt = dx_dt';

function [dx_dt]= eqn1(t,x)
%a function which returns a rate of change vector
% global AlpaIme1 AlpaIme2 AlpaCln3 AlpaRpd3Sin3 BetaIme1 BetaIme2 BetaCln3 n1 n2 n3 n4 n5 n6 n7 n8 KIme11 KIme22 KN2 KIme12 KRpd3 GammaPIme2 KPIme2 dIme2 dIme1 KN2Cln3 GammaPCdk1 KPIme2Cdk1 dCln3 dRpd3  N2 KCln3Ime1 n9 GammaPIme1Cln3 n10 KPIme1Cln3 KRpd3Ime1 n12 n18 xo yo po zo no ;
  parameters
  dx_dt(1) = (AlpaIme1/xo)+(BetaIme1/xo)*(((x(1))^n1)/((KIme11/xo)^n1+(x(1))^n1))*(KN2^(n2)/(KN2^n2+(N2*no)^n2))*((KCln3Ime1^(n9))/(KCln3Ime1^n9+(x(3)*zo)^n9))*((KRpd3Ime1^(n11))/(KRpd3Ime1^n11+((n12-(n12/(1+exp(t-6))))*po)^n11))-(GammaPIme1Cln3/xo)*(((x(3))^n10)/((KPIme1Cln3/zo)^n10+(x(3))^n10))-(dIme1*x(2)*yo*x(1)) ;

%equation of Ime2 protein
dx_dt(2) = (AlpaIme2/yo)+(BetaIme2/yo)*((((x(2))^n3)/((KIme22/yo)^n3+(x(2))^n3))*(((x(1))^n4)/((KIme12/xo)^n4+(x(1))^n4))*(((KRpd3)^n5)/((KRpd3)^n5+(x(4)*po)^n5)))-(GammaPIme2/yo)*(((x(3))^n6)/((KPIme2/zo)^n6+(x(3))^n6))-dIme2*x(2);

%equation of Cdk1/Cln3 complex
dx_dt(3) = (AlpaCln3/zo)+(BetaCln3/zo)*((N2))^(n7)/((KN2Cln3/no)^n7+(N2)^n7)-(GammaPCdk1/zo)*(((x(2))^n8)/((KPIme2Cdk1/yo)^n8+(x(2))^n8))-dCln3*x(3);

%equation of Rpd3/Sin3 complex
dx_dt(4) = (AlpaRpd3Sin3/po)-dRpd3*((x(1)*xo)*1/po)*((1/(N2*no))*((n13)/(n15+(exp((1.6*t))))+n14)); 


  
  
 dx_dt = dx_dt';
  
  