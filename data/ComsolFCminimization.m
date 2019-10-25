AngleFinalp=zeros(35,16);
AngleFinal=zeros(35,16);
BxxFinalp=zeros(1,16);
BxxFinal=zeros(1,16);
for nturnsp=10;


n=20; %number of fourier coefficients (2*n+1=number of coefficients)

%%%%%%% Generate zeroth order cos theta coil   %%%%%%%%%%%%%%%%

R=0.3165;%0.3365;%81.86*2.54/2/100;%0.2965;     %appears to be fullscale b0 -> 
L=0.35;              %100*2.54/100;%*8/6;
nturns=nturnsp; %1/4 of the coil (# of half turns)
qt=nturns;
Lx=.40;%0.076;
Ly=.076;%0.40;
%distortions
    beta=0;%-5.2088*1E-6;%0.000005;%0.00005;% %linear distortion (in current shape, not field)
    c=0;%1.6*1E-6;%0.0000006;%;0.00018;  %quadratic distortion (in current shape, not field)
    alpha=R/(qt+1/2)+beta*qt^2/(qt+1/2)/2+c*qt^3/(qt+1/2)/3;
    %delx=R/2/qt:R/qt:R;
   dely=zeros(1,nturns);
   
    %delx=R/2/qt:R/qt:R;
    for i = 1:nturns
    dely(i)=alpha/2+alpha*(i-1)-beta*(i-1)^2/2-c*(i-1)^3/3;
    end
    ang=pi/2+asin(dely/R);
%     xpos=abs(R*cos(ang));
%     ypos=abs(R*sin(ang));
%     R2=(R-2*2.54/100);
%     xend=R2*cos(ang);
%     yend=R2*sin(ang);
    clear str;
    clear strp;

%%%%%%%%%%%%%%%%%Freespace Fourier Coil%%%%%%%%%%%%%%%%%%%%%%5    
    
%%% This is fitter garbage. 
%options = optimoptions('fminsearch','MaxFunEvals',100000,'TolFun',1E-10,'TolX',1E-40);
options = optimset('MaxIter',10000,'MaxFunEvals',10000,'TolFun',1E-20,'TolX',1E-20);
%maxwell fit error function changed to a single variable 'coeff'.

minimize=@(angle)PositionField(angle,R,Lx,Ly,n);
x0=ang;
%minimize=@(coeff)mwfiterr(coeff,coords,coords,coords,bx,by,bz,bxfit.ModelTerms,mwcoefftransform);



[AngleFinalpp,BxxFinalpp,Bxxexitflagpp,BxxGOFpp]=fminunc(minimize,x0,options);
%%%END OF Freespace Fourier COil
%%% begin minimization of Fourier coil with Boundaries, fminunc with comsol


%MiGrad minimization
options = optimset('MaxIter',10000,'MaxFunEvals',100000,'TolFun',1E-20,'TolX',1E-20);
minimize=@(angle)ComsolFourierCoil(angle,[0,0]);
x0=AngleFinalpp;%(1:33,19)';
[AngleFinal,BxxFinal,Bxxexitflag,BxxGOF]=fminunc(minimize,x0,options);
%disp(['Finished Migrad for ', num2str(nturnsp),', with Bxx=' num2str(BxxFinal(nturnsp)),', starting detailed search'])
end


% %simplex minimization (time consuming but better)
% x0=AngleFinalp(1:nturns,nturnsp);
% options = optimset('MaxIter',1000,'MaxFunEvals',1000,'TolFun',1E-20,'TolX',1E-20);
% [AngleFinal(1:nturns,nturnsp),BxxFinal(nturnsp),Bxxexitflag,BxxGOF]=fminsearch(minimize,x0,options);
% 
% disp(['Finished 1000 step search for ',num2str(nturnsp),' # of coils']);

%end