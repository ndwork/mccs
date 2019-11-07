
function [recon,senseMaps] = mri_iSENSEnn( kData, lambda_x, lambda_s, varargin )
  % recon = mri_iSENSEnn( kData, lambda [, 'maxOuterIter', maxOuterIter ] )
  %
  % Inputs:
  % kData is an array of size ( Ny, Nx, nSlices, nCoils ) of kSpace values
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
  p.addParameter( 'maxOuterIter', 3, @ispositive );
  p.addParameter( 'outDir', './out', @(x) true );
  p.addParameter( 'range', [], @isnumeric );
  p.addParameter( 'verbose', false, @islogical );
  p.parse( kData, lambda_x, lambda_s, varargin{:} );
  doCheckAdjoint = p.Results.doCheckAdjoint;
  maxOuterIter = p.Results.maxOuterIter;
  outDir = p.Results.outDir;
  range = p.Results.range;
  verbose = p.Results.verbose;

  ssqRecon = mri_ssqRecon( kData );  % (Ny, Nx, nSlices, nCoils )
  recon = ssqRecon;

  if verbose == true && ~exist( outDir, 'dir' ), mkdir( outDir ); end

  roughMaps = mccs_makeInitialSensitivityMap( kData );
  senseMaps = roughMaps;

  if verbose == true
    mapsFig = figure;
    senseReconsFig = figure;
    reconFig = figure;
  end

  for iter = 1 : maxOuterIter
    disp([ 'Working on iteration ', num2str(iter), ' of ', num2str(maxOuterIter) ]);

    % Estimate the reconstructed image
    recon = mri_nnReconFISTA_multiCoilMultiSlice( kData, senseMaps, lambda_x, ...
      'nIter', 100, 'initialGuess', recon, 'verbose', verbose );

    if verbose == true
      figure( reconFig ); imshowscale( abs(recon), 5, 'range', range );
      titlenice([ 'recon ', num2str(iter) ]);  drawnow;
      if numel( outDir ) > 0
        saveas( reconFig, [outDir, '/recon_', indx2str(iter,maxOuterIter), '.png'] );
        save( [outDir, '/mat_recon_', indx2str(iter,maxOuterIter), '.mat'], 'recon' );
      end
    end

    % Determine the sensitivity maps
    if iter == maxOuterIter, break; end
    senseMaps = mri_iSENSEnn_makeSensitivityMaps( recon, kData, lambda_s, ...
      'initialGuess', senseMaps, 'verbose', verbose, 'checkAdjoints', doCheckAdjoint, ...
      'maxIterOpt', 100 );

    if verbose == true
      figure( mapsFig );  showImageCube( abs( senseMaps ), 5 );
      if numel( outDir ) > 0
        saveas( mapsFig, [outDir, '/senseMaps_', indx2str(iter,maxOuterIter), '.png'] );
      end

      senseRecons = bsxfun( @times, senseMaps, recon );
      figure( senseReconsFig );  showImageCube( abs(senseRecons), 5 );
      if numel( outDir ) > 0
        saveas( senseReconsFig, [outDir, '/senseRecons_', indx2str(iter,maxOuterIter), '.png'] );
      end
    end
  end

end









