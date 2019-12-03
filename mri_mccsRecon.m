
function [recon,senseMaps] = mri_mccsRecon( kData, lambda_x, lambda_s, lambda_h, varargin )
  % recon = mri_mccsRecon( kData, lambda [, 'maxOuterIter', maxOuterIter, 'noiseCov', noiseCov ] )
  %
  % Inputs:
  % kData is an array of size ( Ny, Nx, nSlices, nCoils ) of kSpace values
  % lambda_x - regularization parameter for compressed sensing reconstruction
  % lambda_s - regularization parameter for determining sensitivity maps
  %
  % Optional Inputs:
  % maxOuterIter - scalar specifying the number of outer iterations
  % noiseCov - a 2D array of size nCoils x nCoils specifying the noise covariance
  %
  % Outputs:
  % recon - the final reconstructed volue
  %
  % Written by Nicholas Dwork, Copyright 2019
  %
  % https://github.com/ndwork/dworkLib.git
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  p = inputParser;
  p.addRequired( 'kData', @isnumeric );
  p.addRequired( 'lambda_x', @(x) numel(x) == 1 && x >= 0 );
  p.addRequired( 'lambda_s', @(x) numel(x) == 1 && x >= 0 );
  p.addParameter( 'doCheckAdjoint', false, @islogical );
  p.addParameter( 'kcf', 0.0035, @(x) x >= 0 && numel(x) == 1 );
  p.addParameter( 'maxOuterIter', 30, @ispositive );
  p.addParameter( 'noiseCov', [], @isnumeric );
  p.addParameter( 'outDir', './out', @(x) true );
  p.addParameter( 'range', [], @isnumeric );
  p.addParameter( 'verbose', false, @islogical );
  p.parse( kData, lambda_x, lambda_s, varargin{:} );
  doCheckAdjoint = p.Results.doCheckAdjoint;
  kcf = p.Results.kcf;
  maxOuterIter = p.Results.maxOuterIter;
  noiseCov = p.Results.noiseCov;
  outDir = p.Results.outDir;
  range = p.Results.range;
  verbose = p.Results.verbose;

  ssqRecon = mri_ssqRecon( kData );  % (Ny, Nx, nSlices, nCoils )
  recon = ssqRecon;
  %sakeEspiritL1Recon = sakeRecon2D( kData, 'type', 'espiritL1' );  % (Ny, Nx, nSlices, nCoils )
  %recon = sakeEspiritL1Recon;

  [ Ny, Nx ] = size( recon );
  nPix = Ny * Nx;

  if verbose == true && ~exist( outDir, 'dir' ), mkdir( outDir ); end

  roughMaps = mccs_makeInitialSensitivityMap( kData );
  %roughMaps = mri_makeSensitivityMaps( kData );
  senseMaps = roughMaps;

  if verbose == true
    mapsFig = figure;
    senseReconsFig = figure;
    reconFig = figure;
    logID = fopen( [outDir, '/mccs.log' ], 'w' );
    fprintf( logID, 'mdm, mdm-2, niqe, piqe, maxSenseMapMap\n' );
  end

  % Check to see if there was previous processing to take advantage of
  if exist( [outDir, '/mat_recon.mat'], 'file' )
    load( [outDir, '/mat_senseMaps.mat'], 'senseMaps', 'iter' );
    senseIter = iter;
    load( [outDir, '/mat_recon.mat'], 'recon', 'iter' );
    if iter ~= senseIter, error( 'Something went wrong with saving' ); end
  else
    iter = 1;
  end

  while iter < maxOuterIter
    disp([ 'Working on iteration ', num2str(iter), ' of ', num2str(maxOuterIter) ]);

    % Determine the sensitivity maps
    senseMaps = mccs_makeSensitivityMaps( recon, kData, kcf, lambda_s, lambda_h, ...
      'initialGuess', senseMaps, 'noiseCov', noiseCov, 'maxIterOpt', 90, ...
      'verbose', verbose, 'doCheckAdjoint', doCheckAdjoint );

    if verbose == true
      figure( mapsFig );  showImageCube( abs( senseMaps ), 5 );
      if numel( outDir ) > 0
        saveas( mapsFig, [outDir, '/senseMaps_', indx2str(iter,maxOuterIter), '.png'] );
        if mod( iter, 10 ) == 0
          save( [outDir, '/mat_senseMaps_', indx2str(iter,maxOuterIter), '.mat'], 'senseMaps' );
        end
      end

      senseRecons = bsxfun( @times, senseMaps, recon );
      figure( senseReconsFig );  showImageCube( abs(senseRecons), 5 );
      if numel( outDir ) > 0
        saveas( senseReconsFig, [outDir, '/senseRecons_', indx2str(iter,maxOuterIter), '.png'] );
      end
    end

    % Estimate the reconstructed image
    recon = mri_csReconFISTA_multiCoilMultiSlice( kData, senseMaps, lambda_x*nPix, ...
      'nIter', 30, 'initialGuess', recon, 'noiseCov', noiseCov, 'verbose', verbose );

    maxAbsRecon = max( abs( recon(:) ) );
    if verbose == true
      figure( reconFig ); imshowscale( abs(recon), 5, 'range', range );
      titlenice([ 'recon ', num2str(iter) ]);  drawnow;
      if numel( outDir ) > 0
        saveas( reconFig, [outDir, '/recon_', indx2str(iter,maxOuterIter), '.png'] );
        save( [outDir, '/mat_senseMaps.mat'], 'senseMaps', 'iter' );
        save( [outDir, '/mat_recon.mat'], 'recon', 'iter' );
        if mod( iter, 10 ) == 0
          save( [outDir, '/mat_recon_', indx2str(iter,maxOuterIter), '.mat'], 'recon' );
        end
      end

      % Calculate image statistics
      maxSenseMapMag = max( abs( senseMaps(:) ) );
      
      % Calculate quality scores
      if maxAbsRecon == 0
        normAbsRecon = recon;
      else
        normAbsRecon = abs( recon ) / max( abs( recon(:) ) );
      end
      mdmScore = calcMdmMetric( normAbsRecon );
      mdmScore2 = calcMdmMetric( 1-normAbsRecon );
      niqeScore = niqe( normAbsRecon );
      piqeScore = piqe( normAbsRecon );
      fprintf( logID, '%f, %f, %f, %f, %f\n', mdmScore, mdmScore2, ...
        niqeScore, piqeScore, maxSenseMapMag );
    end
    if maxAbsRecon == 0, return; end

    iter = iter + 1;
  end

end


