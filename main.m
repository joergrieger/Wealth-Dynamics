clear
global exstate endstate1 xgrid
global FL FH GL GH
global gamma1 gamma2 beta
global div1 end1 end2
global passet11 passet12 passet13 passet14 passet15 passet16 passet17 passet18
global ch11 ch12 ch13 ch14 ch15 ch16 ch17 ch18
global ch21 ch22 ch23 ch24 ch25 ch26 ch27 ch28
global Value101 Value102 Value103 Value104 Value105 Value106 Value107 Value108
global Value201 Value202 Value203 Value204 Value205 Value206 Value207 Value208
global ev1 ev2 psi1 psi2 growth
global m
nendog=101;
%---------------------------------------------------------------------------
% Variables for the time iteration
%---------------------------------------------------------------------------
passet1n    =ones(nendog,8)*(-200);
pbondn      =zeros(nendog,8);
hh1asset1n  =zeros(nendog,8);
hh1bondn    =zeros(nendog,8);
ch1n        =zeros(nendog,8);
ch2n        =zeros(nendog,8);
Value1n     =ones(nendog,8);
Value2n     =ones(nendog,8);

passet1     =ones(nendog,8)*(-200);
pbond       =zeros(nendog,8);
hh1asset1   =zeros(nendog,8);
hh1bond     =zeros(nendog,8);
ch1         =ones(nendog,8);
ch2         =ones(nendog,8);
Value1      =ones(nendog,8)*0;
Value2      =ones(nendog,8)*0;
optimSolutions=ones(nendog,16,18)*(-1.0);
%---------------------------------------------------------------------------
% Parameterize the Economy
%---------------------------------------------------------------------------

%----------------------------------
% Beliefs
%----------------------------------

a1=0.50;
a2=0.14;
a3=0.14;
a4=0.14;

phix=0.43;
alpha1=0.57;
alpha2=0.57;
etah=1.60;
eta2=1.00;

A=[a1, alpha1-a1, alpha2-a1, 1+a1-alpha1-alpha2;
   a2, alpha1-a2, alpha2-a2, 1+a2-alpha1-alpha2;
   a3, alpha1-a3, alpha2-a3, 1+a3-alpha1-alpha2;
   a4, alpha1-a4, alpha2-a4, 1+a4-alpha1-alpha2];

transMat=[0.43 0.57;
          0.57 0.43];

Gamma=[transMat(1,1)*A transMat(1,2)*A;
       transMat(2,1)*A transMat(2,2)*A];


FH=[   etah*transMat(1,1)*A transMat(1,2)*(1-etah*(transMat(1,1)))/(transMat(1,2))*A;
       etah*transMat(2,1)*A transMat(2,2)*(1-etah*(transMat(2,1)))/(transMat(2,2))*A];

FL=1/(1-alpha1)*(Gamma-alpha1*FH);
GH=[   eta2*transMat(1,1)*A transMat(1,2)*(1-eta2*(transMat(1,1)))/(transMat(1,2))*A;
       eta2*transMat(2,1)*A transMat(2,2)*(1-eta2*(transMat(2,1)))/(transMat(2,2))*A];

GL=1/(1-alpha1)*(Gamma-alpha1*GH);

%--------------------------------------------
% Preferences
%-------------------------------------------
gamma1=1.5;
gamma2=1.5;
psi1=1/1.5;
psi2=1/1.5;
beta=0.96;
Value1=log((1-beta)*0.5)*Value1;
Value2=log((1-beta)*0.5)*Value2;
%-----------------------------------------------
% Other primitives of the Economy
%-----------------------------------------------
growth=[1.054 1.054 1.054 1.054  0.982 0.982 0.982 0.982];
div1=[0.15 0.15 0.15 0.15 0.15 0.15 0.15 0.15];
m1=1.00;
m2=1.00;
m=[m1 m1 m1 m1 m2 m2 m2 m2]*1.00; % margin requirements
eh=0.50;
el=1-eh;
end1=[eh eh eh eh el el el el].*(1-div1);
xgrid=linspace(0.00,1.000,nendog);
%--------------------------------------------------
% Set parameters for solvers
%--------------------------------------------------
opts=optiset('solver','ipopt','maxiter',10000,'display','iter','tolafun',1e-7,'tolrfun',1e-7);
options = optimoptions(@fmincon,'MaxIter',50000,'MaxFunEvals',5000000,'Display','off','Algorithm','interior-point','SubproblemAlgorithm','ldl-factorization','FinDiffType','forward','TolCon',1e-8,'TolX',1e-20);

%----------------------------------------------
% Load previous Solutions
%-----------------------------------------------
%load('Solutions_diff_3_m100_p150_g015_e160.mat');
load('Solutions_tempa.mat');
%---------------------------------------
% Start time iteration
%---------------------------------------
max_error=100;
time_iter=0;

while max_error>1e-4
    time_iter=time_iter+1;
    disp(time_iter);
    passet11=pchip(xgrid,squeeze(passet1(:,1)));
	  passet12=pchip(xgrid,squeeze(passet1(:,2)));
	  passet13=pchip(xgrid,squeeze(passet1(:,3)));
	  passet14=pchip(xgrid,squeeze(passet1(:,4)));
	  passet15=pchip(xgrid,squeeze(passet1(:,5)));
	  passet16=pchip(xgrid,squeeze(passet1(:,6)));
	  passet17=pchip(xgrid,squeeze(passet1(:,7)));
	  passet18=pchip(xgrid,squeeze(passet1(:,8)));

    ch11=pchip(xgrid,squeeze(ch1(:,1)));
	  ch12=pchip(xgrid,squeeze(ch1(:,2)));
	  ch13=pchip(xgrid,squeeze(ch1(:,3)));
	  ch14=pchip(xgrid,squeeze(ch1(:,4)));
	  ch15=pchip(xgrid,squeeze(ch1(:,5)));
	  ch16=pchip(xgrid,squeeze(ch1(:,6)));
	  ch17=pchip(xgrid,squeeze(ch1(:,7)));
    ch18=pchip(xgrid,squeeze(ch1(:,8)));

    ch21=pchip(xgrid,squeeze(ch2(:,1)));
	  ch22=pchip(xgrid,squeeze(ch2(:,2)));
	  ch23=pchip(xgrid,squeeze(ch2(:,3)));
	  ch24=pchip(xgrid,squeeze(ch2(:,4)));
    ch25=pchip(xgrid,squeeze(ch2(:,5)));
	  ch26=pchip(xgrid,squeeze(ch2(:,6)));
	  ch27=pchip(xgrid,squeeze(ch2(:,7)));
	  ch28=pchip(xgrid,squeeze(ch2(:,8)));

    Value101=pchip(xgrid,squeeze(Value1(:,1)));
    Value102=pchip(xgrid,squeeze(Value1(:,2)));
    Value103=pchip(xgrid,squeeze(Value1(:,3)));
    Value104=pchip(xgrid,squeeze(Value1(:,4)));
    Value105=pchip(xgrid,squeeze(Value1(:,5)));
    Value106=pchip(xgrid,squeeze(Value1(:,6)));
    Value107=pchip(xgrid,squeeze(Value1(:,7)));
    Value108=pchip(xgrid,squeeze(Value1(:,8)));

    Value201=pchip(xgrid,squeeze(Value2(:,1)));
    Value202=pchip(xgrid,squeeze(Value2(:,2)));
    Value203=pchip(xgrid,squeeze(Value2(:,3)));
    Value204=pchip(xgrid,squeeze(Value2(:,4)));
    Value205=pchip(xgrid,squeeze(Value2(:,5)));
    Value206=pchip(xgrid,squeeze(Value2(:,6)));
    Value207=pchip(xgrid,squeeze(Value2(:,7)));
    Value208=pchip(xgrid,squeeze(Value2(:,8)));

	for exstate=1:8
        disp(exstate);
		    xguess=squeeze(optimSolutions(1,exstate,:));
        %if exstate>1
            %xguess=squeeze(optimSolutions(1,exstate-1,:));
        %end
		    for endstate1=1:nendog
            %if time_iter>1
               xguess=squeeze(optimSolutions(endstate1,exstate,:));
            %end
            if time_iter<3
                disp(endstate1);
            end
            %[x,z,exitflag]=opti_fmincon(@FDummy,xguess,[],[],[],[],[],[],@EulerEquation,opts);
            %xguess=x;
            [x,z,exitflag]=fmincon(@FDummy,xguess,[],[],[],[],[],[],@EulerEquation,options);
            xguess=x;
			      passet1n(endstate1,exstate)=x(17);
			      pbondn(endstate1,exstate)=x(18);
			      hh1asset1n(endstate1,exstate)=x(9);
			      hh1bondn(endstate1,exstate)=x(10);
            ch1n(endstate1,exstate)=x(11);
            ch2n(endstate1,exstate)=x(12);
            Value1n(endstate1,exstate)=log(ev1);
            Value2n(endstate1,exstate)=log(ev2);
            optimSolutions(endstate1,exstate,:)=x;
            %pause;
        end
        %pause;
    end
    diff_asset=max(max(abs(exp(passet1n)-exp(passet1))));
	  diff_bond = max(max(abs(exp(pbondn)-exp(pbond))));
	  diff_hh1a1= max(max(abs((hh1asset1n)-(hh1asset1))));
	  diff_hh1b1= max(max(abs((hh1bondn)-(hh1bond))));
	  diff=[diff_asset diff_bond diff_hh1a1 diff_hh1b1];
	  max_error=max(diff);
	  disp('Maximum Differences in Prices and Portfolios');
	  disp(diff);
    %pause;
	  passet1=passet1n;
	  pbond=pbondn;
	  ch1=ch1n;
	  ch2=ch2n;
    Value1=Value1n;
    Value2=Value2n;
	  hh1asset1=hh1asset1n;
	  hh1bond=hh1bondn;
    save('Solutions_tempa.mat','passet1','pbond','hh1asset1','hh1bond','ch1','ch2','Value1','Value2','optimSolutions','max_error','diff');
    %pause;
    %options = optimset('Display','iter','TolCon',1e-10,'TolX',1e-30,'MaxIter',50000);
end
save('Solutions_diff_2_m100_p066_g015_e160.mat','passet1','pbond','hh1asset1','hh1bond','ch1','ch2','Value1','Value2','optimSolutions','max_error');
