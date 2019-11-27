
function printMetrics( logFile, datacase, sampleFraction, nSamples, recon, imgTitle, outDir, varargin )

  p = inputParser;
  p.addOptional( 'trueRecon', [], @isnumeric );
  p.addParameter( 'senseMaps', [], @isnumeric );
  p.parse( varargin{:} );
  trueRecon = p.Results.trueRecon;
  senseMaps = p.Results.senseMaps;

  if max( abs( recon(:) ) ) == 0
    normAbsRecon = recon;
  else
    normAbsRecon = abs( recon ) / max( abs( recon(:) ) );
  end

  mdmScore = calcMdmMetric( normAbsRecon );
  mdmScore2 = calcMdmMetric( 1-normAbsRecon );
  if max( normAbsRecon(:) ) == 0
    niqeScore = 0;
    piqeScore = 0;
  else
    niqeScore = niqe( normAbsRecon );
    piqeScore = piqe( normAbsRecon );
  end

  if numel( trueRecon ) > 0
    mse = norm( abs(recon(:)) - abs( trueRecon(:) ), 2 ).^2 / numel( recon );
  else
    mse = -1;
  end

  logID = fopen( logFile, 'a' );
  fprintf( logID, ['%d, %f, %d, ', imgTitle, ', %f, %f, %f, %f, %f \n'], ...
    datacase, sampleFraction, nSamples, mdmScore, mdmScore2, niqeScore, piqeScore, mse );
  fclose( logID );

  figH = figure; imshowscale( abs(recon), 5 );  title( imgTitle );
  saveas( figH, [ outDir, '/recon_', imgTitle, '.png' ] );  close( figH );

  if numel( senseMaps ) > 0
    mapsFig = figure;  showImageCube( abs( senseMaps ), 5 );
    saveas( mapsFig, [outDir, '/senseMaps_', imgTitle, '.png'] );  close( mapsFig );

    senseRecons = bsxfun( @times, senseMaps, recon );
    senseReconsFig = figure;  showImageCube( abs(senseRecons), 5 );
    saveas( senseReconsFig, [outDir, '/senseRecons_', imgTitle, '.png'] );
    close( senseReconsFig );

    save( [ outDir, '/recon_', imgTitle, '.mat' ], 'recon', 'senseMaps' );

  else

    save( [ outDir, '/recon_', imgTitle, '.mat' ], 'recon' );
  end
end
