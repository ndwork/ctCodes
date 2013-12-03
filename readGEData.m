function out=readGEData( file )

fid = fopen(file,'rb');
temp=fread(fid,1,'float32');        % 1.0
temp=fread(fid,1,'float32');        % 16.0
SDD=fread(fid,1,'float32');         %source to detector mm
S2I=fread(fid,1,'float32');         %source to iso mm
detcellx=fread(fid,1,'float32');     %detector cell in x mm
view1deg=fread(fid,1,'float32');    %first view angle (degree)
DASsampRate=fread(fid,1,'float32'); %DAS sampling rate (Hz)
rotSpeed=fread(fid,1,'float32');    %Rotation speed(sec)
ViewperRot=fread(fid,1,'float32');  %views per rotation
Nviews=fread(fid,1,'float32');      %total views
Nchan=fread(fid,1,'float32');       %number of chan
Isochan=fread(fid,1,'float32');     %iso channel
Nrows=fread(fid,1,'float32');       %number of rows
detcellz=fread(fid,1,'float32');    %detector cell in z mm
temp=fread(fid,1,'float32');        % 0.0
kv=fread(fid,1,'float32');          % kv
mA= fread(fid,1,'float32');         %mA
tilt=fread(fid,1,'float32');        %tilt

sinoArray = fread(fid,[Nchan, inf],'float32').'/2294.5; % Projection data [array of float]

[totalRows, nchans]=size(sinoArray);

sino=zeros(Nviews,Nchan,Nrows);
indx=1;
for i=1:Nviews
    sino(i,:,:) = reshape(sinoArray(indx:indx+Nrows-1,:)',1,Nchan,Nrows);
    indx=indx+Nrows;
end
sino = flipdim(sino,2);

out=struct;
out.mA=mA;
out.fanNumber = 0;
out.fanIncr  = detcellx/Nchan;   
%out.centerOffset 
out.DID = SDD-S2I;
out.SID = S2I;
out.nViews = Nviews;
out.nDets = Nchan;
out.view1=view1deg;
out.detcellz=detcellz;
out.detcellx=detcellx;

out.sino  = sino;





