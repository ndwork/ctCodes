function out = readSiemensData( file )

  % [float] = 32 bit single precision floating point, little endian, IEEE format
  % [int32] = 32 bit integer, little endian
  fid = fopen(file,'rb');
  mA = fread(fid,1,'float32');        % tube current
  fan_number = fread(fid,1,'int32');  % Fan angle number [int32]
  fan_incr = fread(fid,1,'float32');  % Fan angle increment (2*pi/fan_angle_number) [float]
  center_channel_offset = fread(fid,1,'float32'); % Center channel offset [float]
  DID = fread(fid,1,'float32');       % Detector distance (mm from center of rotation) [float]
  SID = fread(fid,1,'float32');       % Source distance (mm from center of rotation) [float]
  Nviews = fread(fid,1,'int32');      % m (number of tube angles) [int32]
  Ndets = fread(fid,1,'int32');       % n (number of fan angles) [int32]
  sino = fread(fid,[Ndets inf],'float32').'/2294.5; % Projection data [array of float]
  sino = fliplr( sino );
  fclose(fid);

  out = struct;
  out.mA = mA;                        % tube current
  out.fanNumber = fan_number;         % Fan angle number [int32]
  out.fanIncr = fan_incr;             % Fan angle increment (2*pi/fan_angle_number) [float]
  out.centerOffset = center_channel_offset; % Center channel offset [float]
  out.DID = DID;                      % Detector distance (mm from center of rotation) [float]
  out.SID = SID;                      % Source distance (mm from center of rotation) [float]
  out.nViews = Nviews;                % m (number of tube angles) [int32]
  out.nDets = Ndets;                  % n (number of fan angles) [int32]
  out.sino = sino;                    % Projection data [array of float]

end
