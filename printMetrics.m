
function printMetrics( logFile, datacase, sampleFraction, nSamples, recon, imgTitle, outDir, varargin )

  p = inputParser;
  p.addOptional( 'trueRecon', [], @isnumeric );
  p.addParameter( 'senseMaps', [], @isnumeric );
  p.addParameter( 'showScale', 3, @ispositive );
  p.addParameter( 'optLogNames', [] );
  p.addParameter( 'optLogValues', [] );  % must be numeric
  p.parse( varargin{:} );
  trueRecon = p.Results.trueRecon;
  senseMaps = p.Results.senseMaps;
  showScale = p.Results.showScale;
  optLogNames = p.Results.optLogNames;
  optLogValues = p.Results.optLogValues;

  if max( abs( recon(:) ) ) == 0
    normAbsRecon = recon;
  else
    normAbsRecon = abs( recon ) / max( abs( recon(:) ) );
  end

  mdmScore = calcMetricMDM( normAbsRecon );
  mdmScore2 = calcMetricMDM( 1-normAbsRecon );
  if max( normAbsRecon(:) ) == 0
    niqeScore = 0;
    piqeScore = 0;
  else
    niqeScore = niqe( normAbsRecon );
    piqeScore = piqe( normAbsRecon );
  end

  if numel( trueRecon ) > 0
    absRecon = abs( recon );
    absTrueRecon = abs( trueRecon );
    absUnitRecon = absRecon / max( absRecon(:) );
    errImg = absUnitRecon - absTrueRecon;
    errFig = figure;  imshowscale( abs(errImg), showScale, 'range', [0 1] );
    title( [ imgTitle, ' absErr' ] );
    saveas( errFig, [ outDir, '/err_', imgTitle, '.png' ] ); close( errFig );
    mse = norm( abs( errImg(:) ), 2 ).^2 / numel( recon );
    mae = sum( abs( errImg(:) ) ) / numel( recon );  % mean absolute error
    ssimValue = ssim( absUnitRecon, absTrueRecon );
    correlation = dotP( absUnitRecon, absTrueRecon ) / ...
      norm( absUnitRecon(:) ) / norm( absTrueRecon(:) );
    angleErr = acos( correlation );
  else
    mse = -1;
    mae = -1;
    ssimValue = -1;
    angleErr = -1;
  end

  logID = fopen( logFile, 'a' );
  if numel( optLogValues ) > 0
    logValuesString = sprintf( '%f, ' , optLogValues );
    disp([ 'Logging datacase ', num2str(datacase), ' / ', imgTitle, ...
      ' / ', num2str(sampleFraction), ' / ', strjoin( optLogNames, ', ' ), ...
      ', ', logValuesString(1:end-1) ]);
    fprintf( logID, ['%d, %f, %d, ', imgTitle, ', %f, %f, %f, %f, %f, %f, %f, %f, %f ', ...
      repmat( ', %f', [1 numel(optLogValues)] ), '\n' ], ...
      datacase, sampleFraction, nSamples, mdmScore, mdmScore2, niqeScore, piqeScore, ...
      mse, mae, ssimValue, correlation, angleErr, optLogValues(:)' );
  else
    disp([ 'Logging datacase ', num2str(datacase), ' / ', imgTitle, ...
      ' / ', num2str(sampleFraction) ]);
    fprintf( logID, ['%d, %f, %d, ', imgTitle, ', %f, %f, %f, %f, %f, %f, %f, %f, %f \n'], ...
      datacase, sampleFraction, nSamples, mdmScore, mdmScore2, niqeScore, piqeScore, ...
      mse, mae, ssimValue, correlation, angleErr );
  end
  fclose( logID );

  figH = figure; imshowscale( abs(recon), showScale );  title( imgTitle );
  saveas( figH, [ outDir, '/recon_', imgTitle, '.png' ] );  close( figH );

  if numel( senseMaps ) > 0
    mapsFig = figure;  showImageCube( abs( senseMaps ), showScale );
    saveas( mapsFig, [outDir, '/senseMaps_', imgTitle, '.png'] );  close( mapsFig );

    senseRecons = bsxfun( @times, senseMaps, recon );
    senseReconsFig = figure;  showImageCube( abs(senseRecons), showScale );
    saveas( senseReconsFig, [outDir, '/senseRecons_', imgTitle, '.png'] );
    close( senseReconsFig );

    save( [ outDir, '/mat_recon_', imgTitle, '.mat' ], 'recon', 'senseMaps' );

  else

    save( [ outDir, '/mat_recon_', imgTitle, '.mat' ], 'recon' );
  end
end
