
function [recon,senseMaps] = mri_mccsRecon( kData, lambda_x, lambda_s, lambda_h, noiseCoords, ...
  varargin )
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
  p.addRequired( 'lambda_h', @(x) numel(x) == 1 && x >= 0 );
  p.addParameter( 'doCheckAdjoint', false, @islogical );
  p.addParameter( 'kcf', 0.0035, @isnumeric );
  p.addParameter( 'maxOuterIter', 30, @ispositive );
  p.addParameter( 'nReconIter', 30, @ispositive );
  p.addParameter( 'nSenseMapsIter', 90, @ispositive );
  p.addParameter( 'noiseCov', [], @isnumeric );
  p.addParameter( 'outDir', './out', @(x) true );
  p.addParameter( 'range', [], @isnumeric );
  p.addParameter( 'showScale', 2, @ispositive );
  p.addParameter( 'verbose', false, @islogical );
  p.parse( kData, lambda_x, lambda_s, lambda_h, varargin{:} );
  doCheckAdjoint = p.Results.doCheckAdjoint;
  kcf = p.Results.kcf;
  maxOuterIter = p.Results.maxOuterIter;
  nReconIter = p.Results.nReconIter;
  nSenseMapsIter = p.Results.nSenseMapsIter;
  noiseCov = p.Results.noiseCov;
  outDir = p.Results.outDir;
  range = p.Results.range;
  showScale = p.Results.showScale;
  verbose = p.Results.verbose;

  ssqRecon = mri_ssqRecon( kData, 'multiSlice', true );  % (Ny, Nx, nSlices, nCoils )
  recon = ssqRecon;

  senseMaps = mccs_makeInitialSensitivityMap( kData, noiseCoords, kcf );

  if verbose == true && ~exist( outDir, 'dir' ), mkdir( outDir ); end

  if verbose == true
    if numel( outDir ) > 0
      fftRecons = mri_fftRecon( kData, 'multiSlice', true );
      fftReconsFig = figure;
      showImageCube( abs( squeeze( fftRecons ) ), 3 );
      saveas( fftReconsFig, [outDir, '/fftRecons.jpg'] );
      close( fftReconsFig );
    end

    mapsFig = figure;
    reconMapsFig = figure;
    senseReconsFig = figure;
    reconFig = figure;
    logID = fopen( [outDir, '/mccs.log' ], 'w' );
    fprintf( logID, 'mdm, mdm-2, niqe, piqe, maxSenseMapMag\n' );
  end

  iter = 1;
  [ Ny, Nx, nSlices, nCoils ] = size( kData );
  while iter <= maxOuterIter
    disp([ 'Working on iteration ', num2str(iter), ' of ', num2str(maxOuterIter) ]);

    % Determine the sensitivity maps
    senseMaps = mccs_makeSensitivityMaps( recon, kData, kcf, lambda_s, lambda_h, ...
      'initialGuess', senseMaps, 'noiseCov', noiseCov, 'maxIterOpt', nSenseMapsIter, ...
      'verbose', verbose, 'doCheckAdjoint', doCheckAdjoint );
    reconSenseMaps = cropData( senseMaps, [ Ny Nx nSlices nCoils ] );

    if verbose == true
      figure( mapsFig );  showImageCube( abs( senseMaps ), showScale );
      if numel( outDir ) > 0
        saveas( mapsFig, [outDir, '/senseMaps_', indx2str(iter,maxOuterIter), '.jpg'] );
      end

      figure( reconMapsFig );  showImageCube( abs( reconSenseMaps ), showScale );
      if numel( outDir ) > 0
        saveas( reconMapsFig, [outDir, '/reconSenseMaps_', indx2str(iter,maxOuterIter), '.jpg'] );
      end
    end

    % Estimate the reconstructed image
    recon = mri_csReconFISTA_multiCoilMultiSlice( kData, reconSenseMaps, lambda_x, ...
      'nIter', nReconIter, 'initialGuess', recon, 'noiseCov', noiseCov, 'verbose', verbose );

    maxAbsRecon = max( abs( recon(:) ) );
    if verbose == true
      senseRecons = bsxfun( @times, reconSenseMaps, recon );
      figure( senseReconsFig );  showImageCube( abs(senseRecons), showScale );
      if numel( outDir ) > 0
        saveas( senseReconsFig, [outDir, '/senseRecons_', indx2str(iter,maxOuterIter), '.png'] );
      end

      figure( reconFig ); imshowscale( abs(recon), showScale, 'range', range );
      titlenice([ 'recon ', num2str(iter) ]);  drawnow;
      if numel( outDir ) > 0
        saveas( reconFig, [outDir, '/recon_', indx2str(iter,maxOuterIter), '.png'] );
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

  senseMaps = reconSenseMaps;
end

