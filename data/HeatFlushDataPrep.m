%IV1 to IVC tube example flush. 
r=2.03;
length=100;
j=1;
xtemp=linspace(-r,r,20);
ytemp=xtemp;
coords=[];
for i =1:20
    
    for ii = 1:20
        tempcoords=[xtemp(i)*ones(200,1),ytemp(ii)*ones(200,1),linspace(0,length,200)'];
        coords=[coords;tempcoords];
    end
end
[u,v,w,T] = mphinterp(HFdata,{'u','v','w','Temp'},'coord',coords');
