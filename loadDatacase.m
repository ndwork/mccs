
function [ data, noiseCoords, kcfs, res, lambda_xs, lambda_ss, lambda_hs, lambda_iSENSEnn_xs, ...
  lambda_iSENSEnn_ss, lambda_sparseSENSEs,lambda_nnSENSEs, lambda_espiritL1s, ...
  trueRecon ] = loadDatacase( datacase, varargin )
  % Returns data in hybrid state
  %
  % Index:
  % datacase - integer specifying datacase
  %
  % Ouptuts:
  % data - array of data in hybrid state
  % res - the intended resolution corresponding to 0.5 in k-space
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

  lambda_iSENSEnn_xs = 1d-3;
  lambda_iSENSEnn_ss = 1d-2;
  lambda_sparseSENSEs = lambda_iSENSEnn_xs;
  lambda_nnSENSEs = lambda_iSENSEnn_xs;

  kcfs = 9;  % 9 cycles per meter, determined by running findSenseMapsBW
  res = 1d-3;  %#ok<NASGU>

  lambda_xs = 1d-10;
  lambda_ss = 1d-4;
  lambda_hs = 1d6;
  lambda_espiritL1s = [ 0.0002 0.001 0.002 0.0025 0.003 0.004 0.02 ];

  trueRecon = [];

  switch datacase

    case 0
      % Simulated datacase
      [ data, trueRecon ] = loadSimData( simCase );

      noiseCoords = [ 1 1 25 25 ];
      res = 0.1d-3;  % meters per pixel
      lambda_xs = [ 1d-8 1d-9 ];
      lambda_ss = [ 1d-4 1d-6 ];
      lambda_hs = [ 1d3 1d7 ];
      lambda_sparseSENSEs = 1d-4;
      lambda_iSENSEnn_xs = 1d-3;
      lambda_iSENSEnn_ss = 1d-5;

    case 1
      if loadData == true
        data = readOldMriDataOrgData( [ dataDir, '/mriData.Org/P14/kspace' ] );
        data = ifft( ifftshift( data, 3 ), [], 3 );
        data = data( :, :, 5:10:end, : );
        data = data( :, :, 4, : );
      end
      noiseCoords = [ 1 1 41 76 ];
      res = 0.5d-3;  % meters per pixel
      lambda_xs = [ 1d-11 1d-10 1d-9 1d-8 1d-7 1d-6 ];
      lambda_ss = [ 1d-6 1d-5 1d-4 1d-3 1d-2 ];
      lambda_hs = [ 1d4 1d5 1d6 1d7 1d8 ];
      lambda_sparseSENSEs = 1d-4;
      lambda_iSENSEnn_xs = 2d-3;  % 1d-3 does almost nothing
      lambda_iSENSEnn_ss = 1d-5;  % 1d-2, 1d-5 is too much;  1d-3, 1d-4 is too much

    case 2
      if loadData == true
        data = readOldMriDataOrgData( [ dataDir, '/mriData.Org/P17/kspace' ] );
        data = ifft( ifftshift( data, 3 ), [], 3 );
        data = data( :, :, 5:10:end, : );
        data = data( :, :, 4, : );
      end
      noiseCoords = [ 1 1 56 72 ];
      res = 0.5d-3;  % meters per pixel
      lambda_xs = [ 1d-8 1d-10 ];
      lambda_ss = [ 1d-3 1d-6 ];
      lambda_hs = [ 1d2 1d7 ];
      lambda_iSENSEnn_xs = 2d-3;
      lambda_iSENSEnn_ss = 1d-5;  %1d-4 is too high, 1d-5 does almost nothing
      lambda_sparseSENSEs = 1d-4;

    case 3
      if loadData == true
        load( './ESPIRIT/data/brain_8ch.mat', 'DATA' );
        data = permute( DATA, [1 2 4 3] );
      end
      noiseCoords = [ 1 1 35 32 ];
      res = 1.0d-3;  % meters per pixel
      lambda_xs = [ 1d-9 1d-8 ];
      lambda_ss = [ 1d-4 1d-6 ];
      lambda_hs = [ 1d2 1d7 ];
      lambda_iSENSEnn_xs = 1d-7;
      lambda_iSENSEnn_ss = 1d-6;

    case 4
      if loadData == true
        load( './ESPIRIT/data/brain_32ch.mat', 'DATA' );
        data = permute( DATA, [1 2 4 3] );
      end
      noiseCoords = [ 1 1 35 32 ];
      res = 1.0d-3;  % meters per pixel

    case 5
      % FOV 25.6 cm, 256 x 256 (1mm in plane, 1.5mm thick slices)
      if loadData == true
        [data,header] = read_MR_rawdata( [ dataDir, '/cartesian3d/brain3d/P73216.7' ] );   %#ok<ASGLU>
        data = squeeze( data );
        data = ifft( ifftshift( data, 3 ), [], 3 );
        data = data(:,:,64,:);
      end
      noiseCoords = [ 1 1 50 50 ];
      res = 1.0d-3;  % meters per pixel
      lambda_xs = [ 1d-6 1d-8 ];
      lambda_ss = [ 1d-5 1d-6 ];
      lambda_hs = [ 1d2 1d7 ];

    case 6
      % FOV 25.6 cm, 256 x 256 (1mm in plane, 1mm thick slices)
      if loadData == true
        [data,header] = read_MR_rawdata( [ dataDir, '/cartesian3d/brain3d/P73728.7' ] );   %#ok<ASGLU>
        data = squeeze( data );
        data = ifft( ifftshift( data, 3 ), [], 3 );
        data = data(:,:,64,:);
      end
      noiseCoords = [ 1 1 50 50 ];
      res = 1.0d-3;  % meters per pixel
      lambda_xs = [ 1d-6 1d-8 ];
      lambda_ss = [ 1d-5 1d-6 ];
      lambda_hs = [ 1d1 1d7 ];

    case 7
      % FOV 20.5, 256 x 256 (0.8mm in plane, 0.8mm thick slices)
      if loadData == true
        [data,header] = read_MR_rawdata( [ dataDir, '/cartesian3d/brain3d/P74240.7' ] );   %#ok<ASGLU>
        data = squeeze( data );
        data = ifft( ifftshift( data, 3 ), [], 3 );
        data = data(:,:,64,:);
      end
      noiseCoords = [ 1 1 45 45 ];
      res = 0.8d-3;  % meters per pixel
      lambda_xs = [ 1d-6 1d-8 ];
      lambda_ss = [ 1d-6 1d-6 ];
      lambda_hs = [ 1d2 1d7 ];
      lambda_iSENSEnn_xs = 1d-3;
      lambda_iSENSEnn_ss = 1d-2;

    case 8
      % FOV 20 cm, 400 x 400 (0.5mm in plane, 0.5mm thick slices)
      if loadData == true
        [data,header] = read_MR_rawdata( [ dataDir, '/cartesian3d/brain3d/P74752.7' ] );   %#ok<ASGLU>
        data = squeeze( data );
        data = ifft( ifftshift( data, 3 ), [], 3 );
        data = data(:,:,64,:);
      end
      noiseCoords = [ 1 1 70 70 ];
      res = 0.5d-3;  % meters per pixel
      lambda_xs = [ 1d-5 1d-8 ];
      lambda_ss = [ 1d-4 1d-6 ];
      lambda_hs = [ 1d2 1d7 ];
      lambda_iSENSEnn_xs = 1d-4;
      lambda_iSENSEnn_ss = 1d-3;
      lambda_nnSENSEs = 1d-2;

    case 9
      % Ankle, 256 x 256
      if loadData == true
        [data,header] = read_MR_rawdata( [ dataDir, '/cartesian3d/body3d/P26112.7' ] );   %#ok<ASGLU>
        data = squeeze( data );
        data = ifft( ifftshift( data, 3 ), [], 3 );
        data = data(:,:,80,:);
      end
      noiseCoords = [ 10 10 80 80 ];
      res = 1.0d-3;  % meters per pixel
      lambda_xs = [ 1d-5 1d-8 ];
      lambda_ss = [ 1d-3 1d-6 ];
      lambda_hs = [ 1d2 1d7 ];
      lambda_iSENSEnn_xs = 1d-7;
      lambda_iSENSEnn_ss = 1d-6;
      lambda_nnSENSEs = 1d-2;

    case 10
      % shoulder, 
      if loadData == true
        [data,header] = read_MR_rawdata( [ dataDir, '/cartesian3d/body3d/P23040.7' ] );   %#ok<ASGLU>
        data = squeeze( data );
        data = ifft( ifftshift( data, 3 ), [], 3 );
        data = data(:,:,245,:);
      end
      noiseCoords = [ 10 10 60 60 ];
      res = 1.0d-3;  % meters per pixel
      lambda_xs = [ 1d-5 1d-8 ];
      lambda_ss = [ 1d-4 1d-6 ];
      lambda_hs = [ 1d2 1d7 ];
      lambda_iSENSEnn_xs = 1d-7;
      lambda_iSENSEnn_ss = 1d-6;
      lambda_nnSENSEs = 1d-3;

    case 11
      % Ankle, 256 x 256
      if loadData == true
        [data,header] = read_MR_rawdata( [ dataDir, '/cartesian3d/body3d/P27648.7' ] );   %#ok<ASGLU>
        data = squeeze( data );
        data = ifft( ifftshift( data, 3 ), [], 3 );
        data = data(:,:,50,:);
      end
      noiseCoords = [ 10 10 100 100 ];
      res = 1.0d-3;  % meters per pixel
      lambda_xs = [ 1d-6 1d-8 ];
      lambda_ss = [ 1d-5 1d-6 ];
      lambda_hs = [ 1d7 1d2 ];
      lambda_iSENSEnn_xs = 1d-7;
      lambda_iSENSEnn_ss = 1d-5;
      lambda_nnSENSEs = 1d-2;

    otherwise
      error( 'This datacase doesn''t exist' );
  end

end
