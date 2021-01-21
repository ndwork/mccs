
function [ data, noiseCoords, kcf, lambda_xs, lambda_ss, lambda_hs, lambda_iSENSEnn_x, ...
  lambda_iSENSEnn_s, lambda_csSENSE,lambda_nnSENSE, trueRecon ] = loadDatacase( ...
    datacase, sampleFractions, varargin )
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

  p = inputParser;
  p.addOptional( 'loadData', true, @islogical );
  p.addParameter( 'simCase', 1, @isnumeric );
  p.parse( varargin{:} );
  loadData = p.Results.loadData;
  simCase = p.Results.simCase;
  data = [];

  dataDir = '/Volumes/NDWORK128GB/';
  if ~exist( dataDir, 'dir' )
    dataDir = '/Users/nicholasdwork/.DellEMS/Vols/disk2s1/';
    if ~exist( dataDir, 'dir' )
      dataDir = '/Users/ndwork/Documents/Data/';
    end
  end

  lambda_iSENSEnn_x = 1d-3;
  lambda_iSENSEnn_s = 1d-2;
  lambda_csSENSE = lambda_iSENSEnn_x;
  lambda_nnSENSE = lambda_iSENSEnn_x;

  kcf = 0.005;
  lambda_xs = 1d-10 * ones( numel( sampleFractions ), 1 );
  lambda_ss = 1d-4 * ones( numel( sampleFractions ), 1 );
  lambda_hs = 1d2 * ones( numel( sampleFractions ), 1 );

  trueRecon = [];

  switch datacase

    case 0
      % Simulated datacase
      kcf = 0.002;

      [ data, trueRecon ] = loadSimData( simCase );

      noiseCoords = [ 1 1 25 25 ];
      lambda_xs = 1d-6 * ones( numel( sampleFractions ), 1 );
      lambda_ss = 1d-1 * ones( numel( sampleFractions ), 1 );
      lambda_hs = 1d2 * ones( numel( sampleFractions ), 1 );
      lambda_csSENSE = 1d-4;
      lambda_iSENSEnn_x = 2d-3;
      lambda_iSENSEnn_s = 1d-5;

    case 1
      if loadData == true
        data = readOldMriDataOrgData( [ dataDir, '/mriData.Org/P14/kspace' ] );
        data = ifft( ifftshift( data, 3 ), [], 3 );
        data = data( :, :, 5:10:end, : );
        data = data( :, :, 4, : );
      end
      noiseCoords = [ 1 1 41 76 ];
      kcf = 0.005;
      lambda_xs = 1d-8 * ones( numel( sampleFractions ), 1 );
      lambda_ss = 1d-4 * ones( numel( sampleFractions ), 1 );
      lambda_hs = 1d2 * ones( numel( sampleFractions ), 1 );
      lambda_csSENSE = 1d-4;
      lambda_iSENSEnn_x = 2d-3;  % 1d-3 does almost nothing
      lambda_iSENSEnn_s = 1d-5;  % 1d-2, 1d-5 is too much;  1d-3, 1d-4 is too much

    case 2
      if loadData == true
        data = readOldMriDataOrgData( [ dataDir, '/mriData.Org/P17/kspace' ] );
        data = ifft( ifftshift( data, 3 ), [], 3 );
        data = data( :, :, 5:10:end, : );
        data = data( :, :, 4, : );
      end
      noiseCoords = [ 1 1 56 72 ];
      kcf = 0.005;
      lambda_xs = 1d-10 * ones( numel( sampleFractions ), 1 );
      lambda_ss = 1d-3 * ones( numel( sampleFractions ), 1 );
      lambda_hs = 1d2 * ones( numel( sampleFractions ), 1 );
      lambda_iSENSEnn_x = 2d-3;
      lambda_iSENSEnn_s = 1d-5;  %1d-4 is too high, 1d-5 does almost nothing
      lambda_csSENSE = 1d-4;

    case 3
      if loadData == true
        load( './ESPIRIT/data/brain_8ch.mat', 'DATA' );
        data = permute( DATA, [1 2 4 3] );
      end
      noiseCoords = [ 1 1 35 32 ];
      kcf = 0.006;
      lambda_xs = 1d-9 * ones( numel( sampleFractions ), 1 );
      lambda_ss = 1d-4 * ones( numel( sampleFractions ), 1 );
      lambda_hs = 1d2 * ones( numel( sampleFractions ), 1 );
      lambda_iSENSEnn_x = 1d-7;
      lambda_iSENSEnn_s = 1d-6;

    case 4
      if loadData == true
        load( './ESPIRIT/data/brain_32ch.mat', 'DATA' );
        data = permute( DATA, [1 2 4 3] );
      end
      noiseCoords = [ 1 1 35 32 ];

    case 5
      % FOV 25.6 cm, 256 x 256 (1mm in plane, 1.5mm thick slices)
      if loadData == true
        [data,header] = read_MR_rawdata( [ dataDir, '/cartesian3d/brain3d/P73216.7' ] );   %#ok<ASGLU>
        data = squeeze( data );
        data = ifft( ifftshift( data, 3 ), [], 3 );
        data = data(:,:,64,:);
      end
      noiseCoords = [ 1 1 50 50 ];
      kcf = 0.003;
      lambda_xs = 1d-6 * ones( numel( sampleFractions ), 1 );
      lambda_ss = 1d-5 * ones( numel( sampleFractions ), 1 );
      lambda_hs = 1d2 * ones( numel( sampleFractions ), 1 );

    case 6
      % FOV 25.6 cm, 256 x 256 (1mm in plane, 1mm thick slices)
      if loadData == true
        [data,header] = read_MR_rawdata( [ dataDir, '/cartesian3d/brain3d/P73728.7' ] );   %#ok<ASGLU>
        data = squeeze( data );
        data = ifft( ifftshift( data, 3 ), [], 3 );
        data = data(:,:,64,:);
      end
      noiseCoords = [ 1 1 50 50 ];
      kcf = 0.004;
      lambda_xs = 1d-6 * ones( numel( sampleFractions ), 1 );
      lambda_ss = 1d-5 * ones( numel( sampleFractions ), 1 );
      lambda_hs = 1d1 * ones( numel( sampleFractions ), 1 );

    case 7
      % FOV 20.5, 256 x 256 (0.8mm in plane, 0.8mm thick slices)
      if loadData == true
        [data,header] = read_MR_rawdata( [ dataDir, '/cartesian3d/brain3d/P74240.7' ] );   %#ok<ASGLU>
        data = squeeze( data );
        data = ifft( ifftshift( data, 3 ), [], 3 );
        data = data(:,:,64,:);
      end
      noiseCoords = [ 1 1 45 45 ];
      kcf = 0.008;
      lambda_xs = 1d-6 * ones( numel( sampleFractions ), 1 );
      lambda_ss = 1d-6 * ones( numel( sampleFractions ), 1 );
      lambda_hs = 1d2 * ones( numel( sampleFractions ), 1 );
      lambda_iSENSEnn_x = 1d-3;
      lambda_iSENSEnn_s = 1d-2;

    case 8
      % FOV 20 cm, 400 x 400 (0.5mm in plane, 0.5mm thick slices)
      if loadData == true
        [data,header] = read_MR_rawdata( [ dataDir, '/cartesian3d/brain3d/P74752.7' ] );   %#ok<ASGLU>
        data = squeeze( data );
        data = ifft( ifftshift( data, 3 ), [], 3 );
        data = data(:,:,64,:);
      end
      noiseCoords = [ 1 1 70 70 ];
      kcf = 0.005;
      lambda_xs = 1d-5 * ones( numel( sampleFractions ), 1 );
      lambda_ss = 1d-4 * ones( numel( sampleFractions ), 1 );
      lambda_hs = 1d2 * ones( numel( sampleFractions ), 1 );
      lambda_iSENSEnn_x = 1d-4;
      lambda_iSENSEnn_s = 1d-3;
      lambda_nnSENSE = 1d-2;

    case 9
      % Ankle, 256 x 256
      if loadData == true
        [data,header] = read_MR_rawdata( [ dataDir, '/cartesian3d/body3d/P26112.7' ] );   %#ok<ASGLU>
        data = squeeze( data );
        data = ifft( ifftshift( data, 3 ), [], 3 );
        data = data(:,:,80,:);
      end
      noiseCoords = [ 10 10 80 80 ];
      kcf = 0.005;
      lambda_xs = 1d-5 * ones( numel( sampleFractions ), 1 );
      lambda_ss = 1d-3 * ones( numel( sampleFractions ), 1 );
      lambda_hs = 1d2 * ones( numel( sampleFractions ), 1 );
      lambda_iSENSEnn_x = 1d-7;
      lambda_iSENSEnn_s = 1d-6;
      lambda_nnSENSE = 1d-2;

    case 10
      % shoulder, 
      if loadData == true
        [data,header] = read_MR_rawdata( [ dataDir, '/cartesian3d/body3d/P23040.7' ] );   %#ok<ASGLU>
        data = squeeze( data );
        data = ifft( ifftshift( data, 3 ), [], 3 );
        data = data(:,:,245,:);
      end
      noiseCoords = [ 10 10 60 60 ];
      kcf = 0.005;
      lambda_xs = 1d-5 * ones( numel( sampleFractions ), 1 );
      lambda_ss = 1d-4 * ones( numel( sampleFractions ), 1 );
      lambda_hs = 1d2 * ones( numel( sampleFractions ), 1 );
      lambda_iSENSEnn_x = 1d-7;
      lambda_iSENSEnn_s = 1d-6;
      lambda_nnSENSE = 1d-3;

    case 11
      % Ankle, 256 x 256
      if loadData == true
        [data,header] = read_MR_rawdata( [ dataDir, '/cartesian3d/body3d/P27648.7' ] );   %#ok<ASGLU>
        data = squeeze( data );
        data = ifft( ifftshift( data, 3 ), [], 3 );
        data = data(:,:,50,:);
      end
      noiseCoords = [ 10 10 100 100 ];
      kcf = 0.005;
      lambda_xs = 1d-6 * ones( numel( sampleFractions ), 1 );
      lambda_ss = 1d-5 * ones( numel( sampleFractions ), 1 );
      lambda_hs = 1d2 * ones( numel( sampleFractions ), 1 );
      lambda_iSENSEnn_x = 1d-7;
      lambda_iSENSEnn_s = 1d-5;
      lambda_nnSENSE = 1d-2;

    otherwise
      error( 'This datacase doesn''t exist' );
  end

end
