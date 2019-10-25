
scrsz = get(groot,'ScreenSize');
fig=figure('Position',[1 scrsz(4)/2 scrsz(3)/2 scrsz(4)/1.3]);
%Z = peaks;
%histogram(z(1,:),[0:0.02:1.08]);
%line(x(:,1),y(:,1),z(:,1));


%axis tight manual

derp=false;
start=true;
jj=1;
% 
loops =30; %size(z,1);
clear F2; 
F2(loops) = struct('cdata',[],'colormap',[]);
%  myVideo=VideoWriter('heatflushmovie.avi');
%  myVideo.FrameRate=2.5;
%  open(myVideo);
%  writeVideo(myVideo,F2);
%  close(myVideo);

%windows powerpoint usually freezes and skips half the movie
%when it tries to play a movie
%so I pad the start with ones. 
for j = [1,1,1,1,1,1:loops]
    
    if j>1
        start=false;
    end
    
    if derp == true
        jj=jj+1;
    end
    if j ==1
        derp=true;
    end
    
    clf
    subplot(2,1,1)
    axis([-0.4,0.4,0.0,1.1,-0.4,0.4,])
    view(82,28)
    for i = 1:200
    %histogram(z(j,:),[0:0.02:1.08]);
    if start
    line([x(j,i),x(j,i)+0.01*x(j+1,i)],[z(j,i),z(j,i)+0.01*z(j+1,i)],[y(j,i),y(j,i)+0.01*y(j+1,i)]);
    else
    line(x(j-1:j,i),z(j-1:j,i),y(j-1:j,i));
    end
    
    end
    subplot(2,1,2)
 
    histogram(z(j,:),[0:0.02:1.08]);
    drawnow
    if j == 1
    ax = gca;
    ax.NextPlot = 'replaceChildren';
    end
    F2(jj) = getframe(fig);
    
end
% %save to file
%this fucks up a lot, matlab randomly changes the 
%cdata of the plot even though axis are set and not changed.
%the clear F2 may have helped this. 
 myVideo=VideoWriter('heatflushmovie.avi');
 myVideo.FrameRate=2.5;
 open(myVideo);
 writeVideo(myVideo,F2);
 close(myVideo);
