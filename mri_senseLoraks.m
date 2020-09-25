
function recon = mri_senseLoraks( kData, varargin )

  defaultLambda = 1d-3;
  defaultRank = floor( 55 / 256 * size( kData, 1 ) );

  p = inputParser;
  p.addParameter( 'lambda', defaultLambda, @ispositive );
  p.addParameter( 'rank', defaultRank, @ispositive );
  p.addParameter( 'senseMaps', [], @isnumeric );
  p.addParameter( 'verbose', false, @(x) islogical(x) || isnumeric(x) );
  p.parse( varargin{:} );
  lambda = p.Results.lambda;
  rank = p.Results.rank;
  senseMaps = p.Results.senseMaps;
  verbose = p.Results.verbose;

  if numel( senseMaps ) == 0
    mask = mri_makeIntensityMask( kData );
    senseMaps = mri_makeSensitivityMaps( kData, 'mask', mask, 'verbose', verbose );
  end

  sKData = size( kData );
  nSlices = sKData( 3 );
  recon = zeros( sKData(1:3) );
  for sliceIndx = 1 : nSlices
    sliceKData = squeeze( kData(:,:,sliceIndx,:) );
    sliceSenseMap = squeeze( senseMaps(:,:,sliceIndx,:) );

    kMask = sum( abs( sliceKData ), 3 ) ~= 0;

    sliceRecon = SENSE_LORAKS( sliceKData, kMask, sliceSenseMap, rank, ...
      lambda, [], 'S', 1, [], [] );
    recon(:,:,sliceIndx) = sliceRecon;
  end

end
