function ppB0 = ComsolMSRCoil(pos)
    drool=pos<-1.9;
    if ~isempty(drool)
    for i = 1:length(drool)
    pos(drool)=-1.9+0.01*(i-1);
    end
    end
    drool=pos>1.9;
    if ~isempty(drool)
    for i = 1:length(drool)
    pos(drool)=1.9-0.01*(i-1);
    end
    end
    
    
model=comsolMSR(pos);
    coords=[linspace(-1.8,1.8,100);zeros(1,100);zeros(1,100)];
    
    Bx=mphinterp(model,'mf.Bx','coord',coords);
  
    ppB0=abs(std(Bx)/mean(Bx));
    load('MSRdata.mat');
    if ppB0<ppB0temp
        ppB0temp=ppB0;
        postemp=pos;
        save('MSRdata.mat','postemp','ppB0temp');
    end
    
    disp(ppB0);