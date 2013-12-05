% function MAR
clear
close all

%dataDir = 'C:\Users\ndwork\Desktop\';
dataDir = 'C:\Users\ndwork\Documents\My Stuff\My School\Data\ctMetalArtifact\';
%dataDir = 'C:\Users\Uzair\SkyDrive\Stanford Docs\EE 369C\Project\';
%dataDir = '../data/';
dataFile=[dataDir,'/Siemens_FromEdBoas/precalc_Hep.bin'];

disp(['dataFile: ', dataFile]);

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
    dTheta = 2*dBeta;
    oThetas = [0,dTheta,maxBeta];
    thetas = 0:dTheta:maxBeta-dBeta;
    dSize = 0.001;
    delta = dSize;
    oLines = [-0.3,dSize,0.3];
    sid_m = data.SID / 1000;

    sino = ctFan2Para( data.sino, gammas, betas, sid_m, ...
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
%sizeSino=size(sino);
%[rubbedSino,sinoMask]=rubOutSino(sino,50,thetas,sizeSino(2),dSize, 0, 0, 512, 512, delta, delta, 'Hanning');
%reconRubbed=ctIRadon(rubbedSino, thetas, dSize, 0, 0, 512, 512, delta, delta, 'Hanning');


%%
%recon=ctIRadon(sino, thetas, dSize, 0, 0, 512, 512, delta, delta, 'Hanning');
%figure('name','ctIRadon recon'),imshow(recon,[])

%%
%sizeSino=size(sino);
%reconNonMetal=ctIRadonMetal(sino,thetas,dSize,0,0,512,512,delta,delta,'Hanning');
%figure,imshow(reconNonMetal,[]);

%%
%sizeSino=size(sino);
%marwavRecon=marWav(sino,thetas,sizeSino(2),dSize,0,0,512,512,delta,delta, 'Hanning');
%figure('name','marWav recon'), imshow( marwavRecon, [] );

%%
%sizeSino=size(sino);
%marwavRadRecon=marWavRad(sino,thetas,sizeSino(2),dSize,0,0,512,512,delta,delta, 'Hanning');
%figure('name','marWavRad recon'), imshow( marwavRadRecon, [] );

%%
%sizeSino=size(sino);
%marShearRadRecon=marShearRad(sino,thetas,sizeSino(2),dSize,0,0,512,512,delta,delta, 'Hanning');
%figure('name','marWavRad recon'), imshow( marwavRadRecon, [] );

%%
%ruboutRecon=rubOut(recon,21000);
%figure('name','rubOut recon'), imshow( ruboutRecon, [] );

%%
%reconMDTed=mdt(recon,50);

%%
%sizeSino=size(sino);
%mardsRecon=marDS(sino,thetas,sizeSino(2),dSize,0,0,512,512,delta,delta,'Hanning');

%%
sizeSino=size(sino);
mardsRecon=marPC(sino,thetas,sizeSino(2),dSize,0,0,512,512,delta,delta,'Hanning');

