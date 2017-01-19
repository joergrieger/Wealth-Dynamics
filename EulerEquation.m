function [fvalin,fval] = EulerEquation(x)
% EulerEquation
% Calculates the first order conditions
global exstate endstate1 xgrid
global FL FH GL GH
global gamma1 gamma2 beta
global div1 end1
global passet11 passet12 passet13 passet14 passet15 passet16 passet17 passet18
global ch11 ch12 ch13 ch14 ch15 ch16 ch17 ch18
global ch21 ch22 ch23 ch24 ch25 ch26 ch27 ch28
global Value101 Value102 Value103 Value104 Value105 Value106 Value107 Value108
global Value201 Value202 Value203 Value204 Value205 Value206 Value207 Value208
global ev1 ev2 psi1 psi2 growth
global m 

fval=zeros(18,1);
fvalin=[];
%---------------------------------------------------------
% first set of variables: guesses of future wealth share
%---------------------------------------------------------
hh1wshareguess=zeros(8,1);
for ii=1:8
    hh1wshareguess(ii)=x(ii);
    %hh2wshareguess(ii)=x(16+ii);
end
%-------------------------------------------------------
% guesses for portfolio choices
%-------------------------------------------------------
hh1asset1guess      =x(9);
hh1bondguess        =x(10);
hh2asset1guess      =1-x(9);
hh2bondguess        =-x(10);

%------------------------------------------------------
% guesses for consumption
%------------------------------------------------------
ch1=exp(x(11));
ch2=exp(x(12));
%------------------------------------------------------
% guesses for lagrange multiplier
%------------------------------------------------------
lag1assetshortpos   =max(0,x(13))^2;
lag1assetshortneg   =max(0,-x(13))^2;
lag2assetshortpos   =max(0,x(14))^2;
lag2assetshortneg   =max(0,-x(14))^2;

lag1collpos         =max(0,x(15))^2;
lag1collneg         =max(0,-x(15))^2;
lag2collpos         =max(0,x(16))^2;
lag2collneg         =max(0,-x(16))^2;
%-------------------------------------------------------
% guesses for prices
%------------------------------------------------------
q1g=exp(x(17));
qbg=exp(x(18));
%---------------------------------------------------------
% Interpolate future prices and consumption
%---------------------------------------------------------
ch1t=zeros(8,1);
ch2t=zeros(8,1);
passet1t=zeros(8,1);
Value1t=zeros(8,1);
Value2t=zeros(8,1);

ii=1;
passet1t(ii)    =exp(ppval(passet11,hh1wshareguess(ii)));
ch1t(ii)        =exp(ppval(ch11,hh1wshareguess(ii)));
ch2t(ii)        =exp(ppval(ch21,hh1wshareguess(ii)));
Value1t(ii)     =ppval(Value101,hh1wshareguess(ii));
Value2t(ii)     =ppval(Value201,hh1wshareguess(ii));

ii=2;
passet1t(ii)    =exp(ppval(passet12,hh1wshareguess(ii)));
ch1t(ii)        =exp(ppval(ch12,hh1wshareguess(ii)));
ch2t(ii)        =exp(ppval(ch22,hh1wshareguess(ii)));
Value1t(ii)     =ppval(Value102,hh1wshareguess(ii));
Value2t(ii)     =ppval(Value202,hh1wshareguess(ii));

ii=3;
passet1t(ii)    =exp(ppval(passet13,hh1wshareguess(ii)));
ch1t(ii)        =exp(ppval(ch13,hh1wshareguess(ii)));
ch2t(ii)        =exp(ppval(ch23,hh1wshareguess(ii)));
Value1t(ii)     =ppval(Value103,hh1wshareguess(ii));
Value2t(ii)     =ppval(Value203,hh1wshareguess(ii));

ii=4;
passet1t(ii)    =exp(ppval(passet14,hh1wshareguess(ii)));
ch1t(ii)        =exp(ppval(ch14,hh1wshareguess(ii)));
ch2t(ii)        =exp(ppval(ch24,hh1wshareguess(ii)));
Value1t(ii)     =ppval(Value104,hh1wshareguess(ii));
Value2t(ii)     =ppval(Value204,hh1wshareguess(ii));
ii=5;
passet1t(ii)    =exp(ppval(passet15,hh1wshareguess(ii)));
ch1t(ii)        =exp(ppval(ch15,hh1wshareguess(ii)));
ch2t(ii)        =exp(ppval(ch25,hh1wshareguess(ii)));
Value1t(ii)     =ppval(Value105,hh1wshareguess(ii));
Value2t(ii)     =ppval(Value205,hh1wshareguess(ii));
ii=6;
passet1t(ii)    =exp(ppval(passet16,hh1wshareguess(ii)));
ch1t(ii)        =exp(ppval(ch16,hh1wshareguess(ii)));
ch2t(ii)        =exp(ppval(ch26,hh1wshareguess(ii)));
Value1t(ii)     =ppval(Value106,hh1wshareguess(ii));
Value2t(ii)     =ppval(Value206,hh1wshareguess(ii));

ii=7;
passet1t(ii)    =exp(ppval(passet17,hh1wshareguess(ii)));
ch1t(ii)        =exp(ppval(ch17,hh1wshareguess(ii)));
ch2t(ii)        =exp(ppval(ch27,hh1wshareguess(ii)));
Value1t(ii)     =ppval(Value107,hh1wshareguess(ii));
Value2t(ii)     =ppval(Value207,hh1wshareguess(ii));

ii=8;
passet1t(ii)    =exp(ppval(passet18,hh1wshareguess(ii)));
ch1t(ii)        =exp(ppval(ch18,hh1wshareguess(ii)));
ch2t(ii)        =exp(ppval(ch28,hh1wshareguess(ii)));
Value1t(ii)     =ppval(Value108,hh1wshareguess(ii));
Value2t(ii)     =ppval(Value208,hh1wshareguess(ii));

%--------------------------------------------------------
% First order conditions
%--------------------------------------------------------
Value1t=exp(Value1t);
Value2t=exp(Value2t);
minPA=min(passet1t+div1');
qb1=0;
qb2=0;
q11=0;
q21=0;
if exstate==1 
	transMatrix1=FH;
	transMatrix2=GH;
end
if exstate==2
	transMatrix1=FH;
	transMatrix2=GL;
end
if exstate==3
	transMatrix1=FL;
	transMatrix2=GH;
end
if exstate==4
	transMatrix1=FL;
	transMatrix2=GL;
end
if exstate==5 
	transMatrix1=FH;
	transMatrix2=GH;
end
if exstate==6
	transMatrix1=FH;
	transMatrix2=GL;
end
if exstate==7
	transMatrix1=FL;
	transMatrix2=GH;
end
if exstate==8
	transMatrix1=FL;
	transMatrix2=GL;
end

EValue1=0;
EValue2=0;
for ii=1:8
    qb1 = qb1+transMatrix1(exstate,ii)*Value1t(ii)^(1/psi1-gamma1)*growth(ii)^(-gamma1)*(1-beta)*ch1t(ii)^(-1/psi1);
    qb2 = qb2+transMatrix2(exstate,ii)*Value2t(ii)^(1/psi2-gamma2)*growth(ii)^(-gamma2)*(1-beta)*ch2t(ii)^(-1/psi2);
    
    q11 = q11+transMatrix1(exstate,ii)*Value1t(ii)^(1/psi1-gamma1)*growth(ii)^(1-gamma1)*(1-beta)*ch1t(ii)^(-1/psi1)*(passet1t(ii)+div1(ii));
    q21 = q21+transMatrix2(exstate,ii)*Value2t(ii)^(1/psi2-gamma1)*growth(ii)^(1-gamma2)*(1-beta)*ch2t(ii)^(-1/psi2)*(passet1t(ii)+div1(ii));
    EValue1=EValue1+transMatrix1(exstate,ii)*(Value1t(ii)*growth(ii))^(1-gamma1);
    EValue2=EValue2+transMatrix2(exstate,ii)*(Value2t(ii)*growth(ii))^(1-gamma2);
end
phi1=(1-gamma1)/(1-1/psi1);
phi2=(1-gamma2)/(1-1/psi2);
Value1today=((1-beta)*ch1^(1-1/psi1)+beta*EValue1^(1/phi1))^(phi1/(1-gamma1));
Value2today=((1-beta)*ch2^(1-1/psi2)+beta*EValue2^(1/phi2))^(phi2/(1-gamma2));
MValue1today=(1-beta)*Value1today^(1/psi1)*ch1^(-1/psi1);
MValue2today=(1-beta)*Value2today^(1/psi2)*ch2^(-1/psi2);
ev1=Value1today;
ev2=Value2today;
%--------------------------------------------------------
% Calculate financial wealth share
%-------------------------------------------------------
hh1wsharet=zeros(8,1);
for ii=1:8
    hh1wsharet(ii)=(hh1asset1guess*(passet1t(ii)+div1(ii))+1/growth(ii)*hh1bondguess)/((passet1t(ii)+div1(ii)));
    fval(ii)=hh1wsharet(ii)-hh1wshareguess(ii);
    %fval(16+ii)=1-hh1wshareguess(ii)-hh2wshareguess(ii);
end
%---------------------------------------------------------------------
% First order conditions
%---------------------------------------------------------------------
fval(9)=q1g*MValue1today-(lag1assetshortpos+lag1collpos*(1-m(exstate))*minPA+beta*Value1today^(1/psi1)*EValue1^((1-phi1)/phi1)*q11);
fval(10)=q1g*MValue2today-(lag2assetshortpos+lag2collpos*(1-m(exstate))*minPA+beta*Value2today^(1/psi2)*EValue2^((1-phi2)/phi2)*q21);
fval(11)=qbg*MValue1today-(lag1collpos+beta*Value1today^(1/psi1)*EValue1^((1-phi1)/phi1)*qb1);
fval(12)=qbg*MValue2today-(lag2collpos+beta*Value2today^(1/psi2)*EValue2^((1-phi2)/phi2)*qb2);

%----------------------------------------------------------------------
% Budget constraint
%----------------------------------------------------------------------
income1=end1(exstate)+xgrid(endstate1)*(q1g+div1(exstate))-hh1asset1guess*q1g-hh1bondguess*qbg;
fval(13)=ch1-income1;
fval(14)=1-ch1-ch2;
%-------------------------------------------------------------------
% Complementary slackness conditions
%-------------------------------------------------------------------
fval(15)=lag1assetshortneg-hh1asset1guess;
fval(16)=lag2assetshortneg-hh2asset1guess;
fval(17)=lag1collneg-((1-m(exstate))*(hh1asset1guess*minPA)+hh1bondguess);
fval(18)=lag2collneg-((1-m(exstate))*(hh2asset1guess*minPA)+hh2bondguess);

end

