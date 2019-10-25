%MiGrad minimization
options = optimset('MaxIter',1000,'MaxFunEvals',1000,'TolFun',1E-10,'TolX',1E-10);
minimize=@(position)ComsolMSRCoil(position);
x0=pos;%(1:33,19)';
[PositionFinal,BxxFinal,Bxxexitflag,BxxGOF]=fminsearch(minimize,x0,options);