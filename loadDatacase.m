
function [data,noiseCoords] = loadDatacase( datacase )
  % Returns data in hybrid state
  %
  % Index:
  % datacase - integer specifying datacase
  %
  % Ouptuts:
  % data - array of data in hybrid state
  % noiseCoords - [ minX minY maxX maxY ] specifying noise region in data
  %
  % Written by Nicholas Dwork, Copyright 2019

  dataDir = '/Users/nicholasdwork/.DellEMS/Vols/disk2s1/cartesian3d/';
  if ~exist( dataDir, 'dir' )
    dataDir = '/Users/ndwork/Documents/Data/cartesian3d/';
  end

  switch datacase

    case 1
      data = readOldMriDataOrgData( [ dataDir, '/mriData.Org/P14/kspace' ] );
      data = ifftc( data, [], 3 );
      data = data( :, :, 5:10:end, : );
      data = data( :, :, 4, : );
      noiseCoords = [ 1 1 41 76 ];

    case 2
      data = readOldMriDataOrgData( [ dataDir, '/mriData.Org/P17/kspace' ] );
      data = ifftc( data, [], 3 );
      data = data( :, :, 5:10:end, : );
      data = data( :, :, 4, : );
      noiseCoords = [ 1 1 56 72 ];

    case 3
      load( './ESPIRIT/data/brain_8ch.mat', 'DATA' );
      data = permute( DATA, [1 2 4 3] );
      noiseCoords = [ 1 1 35 32 ];
      %kcf = 0.008;
      %lambda_xs = 1d-12 * ones( numel( sampleFractions ), 1 );
      %lambda_ss = 1d-11 * ones( numel( sampleFractions ), 1 );

    case 4
      load( './ESPIRIT/data/brain_32ch.mat', 'DATA' );
      data = permute( DATA, [1 2 4 3] );
      noiseCoords = [ 1 1 35 32 ];
      %kcf = 0.008;
      %lambda_xs = 1d-12 * ones( numel( sampleFractions ), 1 );
      %lambda_ss = 1d-11 * ones( numel( sampleFractions ), 1 );

    case 5
      % FOV 25.6 cm, 256 x 256 (1mm in plane, 1.5mm thick slices)
      [data,header] = read_MR_rawdata( [ dataDir, '/brain3d/P73216.7' ] );   %#ok<ASGLU>
      data = squeeze( data );
      data = ifftc( data, [], 3 );
      data = data(:,:,64,:);
      noiseCoords = [ 1 1 50 50 ];

    case 6
      % FOV 25.6 cm, 256 x 256 (1mm in plane, 1mm thick slices)
      [data,header] = read_MR_rawdata( [ dataDir, '/brain3d/P73728.7' ] );   %#ok<ASGLU>
      data = squeeze( data );
      data = ifftc( data, [], 3 );
      data = data(:,:,64,:);
      noiseCoords = [ 1 1 50 50 ];

    case 7
      % FOV 20.5, 256 x 256 (0.8mm in plane, 0.8mm thick slices)
      [data,header] = read_MR_rawdata( [ dataDir, '/brain3d/P74240.7' ] );   %#ok<ASGLU>
      data = squeeze( data );
      data = ifftc( data, [], 3 );
      data = data(:,:,64,:);
      noiseCoords = [ 1 1 45 45 ];

    case 8
      % FOV 20 cm, 400 x 400 (0.5mm in plane, 0.5mm thick slices)
      [data,header] = read_MR_rawdata( [ dataDir, '/brain3d/P74752.7' ] );   %#ok<ASGLU>
      data = squeeze( data );
      data = ifftc( data, [], 3 );
      data = data(:,:,64,:);
      noiseCoords = [ 1 1 70 70 ];

    case 9
      % Ankle, 256 x 256
      [data,header] = read_MR_rawdata( [ dataDir, '/body3d/P26112.7' ] );   %#ok<ASGLU>
      data = squeeze( data );
      data = ifftc( data, [], 3 );
      data = data(:,:,80,:);
      noiseCoords = [ 10 10 80 80 ];

    case 10
      % shoulder, 
      [data,header] = read_MR_rawdata( [ dataDir, '/body3d/P23040.7' ] );   %#ok<ASGLU>
      data = squeeze( data );
      data = ifftc( data, [], 3 );
      data = data(:,:,245,:);
      noiseCoords = [ 10 10 60 60 ];

    case 11
      % Ankle, 256 x 256
      [data,header] = read_MR_rawdata( [ dataDir, '/body3d/P27648.7' ] );   %#ok<ASGLU>
      data = squeeze( data );
      data = ifftc( data, [], 3 );
      data = data(:,:,50,:);
      noiseCoords = [ 10 10 100 100 ];

    otherwise
      error( 'This datacase doesn''t exist' );
  end

end
