
function printMetrics( logFile, datacase, sampleFraction, nSamples, recon, imgTitle, varargin )

  p = inputParser;
  p.addOptional( 'trueRecon', [], @isnumeric );
  p.parse( varargin{:} );
  trueRecon = p.Results.trueRecon;

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

end
