
function [recon,senseMaps] = mri_nnSenseRecon( kData, lambda, varargin )
  % [recon,senseMaps] = mri_nnSenseRecon( kData [, 'mask', mask, ...
  %   'senseMaps', senseMaps, 'verbose', true/false ] )
  %
  % Inputs:
  % kData is an array of size ( Ny, Nx, nSlices, nCoils ) of kSpace values
  %
  % https://github.com/ndwork/dworkLib.git
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  p = inputParser;
  p.addRequired( 'kData', @isnumeric );
  p.addParameter( 'senseMaps', [], @isnumeric );
  p.addParameter( 'verbose', false, @islogical );
  p.parse( kData, varargin{:} );
  senseMaps = p.Results.senseMaps;
  verbose = p.Results.verbose;

  if numel( senseMaps ) == 0
    mask = mri_makeIntensityMask( kData );
    senseMaps = mri_makeSensitivityMaps( kData, 'mask', mask, 'verbose', verbose );
  end

  recon = mri_nnReconFISTA_multiCoilMultiSlice( kData, senseMaps, lambda, 'verbose', verbose );

end
