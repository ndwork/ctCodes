% function MAR
clc
clear
close all

dataDir = 'C:\Users\ndwork\Documents\My Stuff\My School\Data\ctMetalArtifact\';
dataFile=[dataDir,'Siemens_FromEdBoas\precalc_Hep.bin'];

sliceIndx=5;

if numel( dataFile ) > 0
    if( regexp(dataFile,'\.bin$') ) %Siemens data
        data = readSiemensData( dataFile );
        
        sSinoF = size( data.sino );
        totFanAngle = data.fanIncr * sSinoF(2);
        gammas = ( data.fanIncr * [1:sSinoF(2)] ) - data.fanIncr*sSinoF(2)/2;
        dBeta = 2*pi/sSinoF(1);
        maxBeta = 2*pi - dBeta;
        betas = 0:dBeta:maxBeta;
        oThetas = [0,dBeta,maxBeta];
        thetas=0:dBeta:maxBeta-dBeta;
        dSize=0.001;
        delta=dSize;
        oLines = [-0.3,dSize,0.3];
        sid_m = data.SID / 1000;
        
        sino = ctFan2para( data.sino, gammas, betas, sid_m, ...
            oThetas, oLines );
    else
        %GE Data
        
        data = readGEData([dataDir,'\',files(i).name]);
        
        sSinoF = size( data.sino );
        totFanAngle = data.fanIncr * sSinoF(2);
        gammas = ( data.fanIncr * [1:sSinoF(2)] ) - data.fanIncr*sSinoF(2)/2;
        dBeta = 2*pi/sSinoF(1);
        maxBeta = 2*pi - dBeta;
        betas = 0:dBeta:maxBeta;
        oThetas = [0,dBeta,maxBeta];
        thetas=0:dBeta:maxBeta-dBeta;
        dSize=0.001;
        delta=dSize;
        oLines = [-0.3,dSize,0.3];
        
        
        sid_m = data.SID / 1000;
        
        sinoF=squeeze(data.sino(:,:,sliceIndx));
        
        sino = ctFan2para( sinoF, gammas, betas, sid_m, ...
            oThetas, oLines );
        
    end
else
    % Simulated Data
    delta=0.001;
    NumDet=500;
    dSize=delta;
    dtheta=1*(pi/180);
    thetas=0:dtheta:pi-dtheta;
    phant=createPhantom(delta, [0, 0], 0.05 );
    sino=ctRadon(phant.im,delta,NumDet,dSize,thetas);
end

%%
recon=ctIRadon(sino, thetas, dSize, 0, 0, 512, 512, delta, delta, 'Hanning');
figure,imshow(recon,[])


%%
% metalMask=findMetal(recon,21000);
sizeSino=size(sino);
reconWav=marWav(sino,thetas,sizeSino(2),dSize,0,0,512,512,delta,delta, 'Hanning');
%%

reconRubbed=rubOut(recon,21000);
%%

reconMDTed=mdt(recon,21000);

