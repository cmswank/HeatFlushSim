

fileID = fopen('/data1/cmswank/SpinSimSwank2/data/spindataTest.dat');

A = fread(fileID, 'double');

%size(A)
Bin=A(1); %timesteps
Event=A(2); %number of spins
nvars=A(3); %number of variables in the datafile. 
A=A(4:end);

B = reshape(A,nvars,Bin,Event);

sx = squeeze(B(1,:,:));

sy = squeeze(B(2,:,:));

sz = squeeze(B(3,:,:));

x =  squeeze(B(4,:,:));

y =  squeeze(B(5,:,:));

z = squeeze(B(6,:,:));

vx =  squeeze(B(7,:,:));

vy =  squeeze(B(8,:,:));

vz = squeeze(B(9,:,:));

bx = squeeze(B(10,:,:));

by = squeeze(B(11,:,:));

bz = squeeze(B(12,:,:));

ex = squeeze(B(13,:,:));

ey = squeeze(B(14,:,:));

ez = squeeze(B(15,:,:));


tlarge = squeeze(B(16,:,:)); 

t = squeeze(tlarge(:,1));

f=(1/t(end):1/t(end):length(t)/t(end))';

Sz=mean(sz,2);
Sx=mean(sx,2);
Sy=mean(sy,2);


%specz=2*fft(mean(sz,2))/length(sz);
%specx=2*fft(mean(sx,2))/length(sx);
%specy=2*fft(mean(sy,2))/length(sy);

%sum(sum(y))
