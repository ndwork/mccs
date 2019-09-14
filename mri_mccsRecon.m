
function [recon,senseMaps] = mri_mccsRecon( kData, lambda, varargin )
  % recon = mri_mccsRecon( kData, lambda [, 'maxOuterIter', maxOuterIter ] )
  %
  % Inputs:
  % kData is an array of size ( Ny, Nx, nSlices, nCoils ) of kSpace values
  %
  % Outputs:
  % recon - the final reconstructed volue
  %
  % https://github.com/ndwork/dworkLib.git
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  p = inputParser;
  p.addRequired( 'kData', @isnumeric );
  p.addRequired( 'lambda', @(x) isnumeric( x ) && numel(x) == 1 );
  p.addParameter( 'maxOuterIter', 4, @ispositive );
  p.addParameter( 'verbose', false, @islogical );
  p.parse( kData, lambda, varargin{:} );
  maxOuterIter = p.Results.maxOuterIter;
  verbose = p.Results.verbose;

  kcf = 0.01;

  ssqRecon = mri_ssqRecon( kData );  % (Ny, Nx, nSlices, nCoils )
  recon = ssqRecon;
  roughMaps = mccs_makeInitialSensitivityMap( kData );
  senseMaps = roughMaps;

  for iter = 1 : maxOuterIter
    disp([ 'Working on iteration ', num2str(iter), ' of ', num2str(maxOuterIter) ]);

    % Determine the sensitivity maps
    senseMaps = mccs_makeSensitivityMaps( recon, kData, kcf, 'initialGuess', senseMaps );

    % Estimate the reconstructed image
    recon = mri_csReconFISTA_multiCoilMultiSlice( kData, senseMaps, lambda, ...
      'verbose', verbose );

  end

  figure; imshowscale( recon, 5 );  titlenice( 'recon' );
  figure; imshowscale( ssqRecon, 5 );  titlenice( 'ssqRecon' );
  figure; showImageCube( senseMaps, 'border', 1, 'borderValue', 'max' );  titlenice( 'sense maps' );
  figure; showImageCube( roughMaps, 'border', 1, 'borderValue', 'max' );  titlenice( 'rough maps' );
end



