
function senseMaps = mccs_makeSensitivityMaps( recon, kData, kcf, varargin )
  % senseMaps = mccs_makeSensitivityMaps( recon, kData, kcf, varargin )
  %
  % Impose a cutoff frequency to restrict the sensitivity maps to low frequencies
  %
  % Inputs:
  % kData - the k-space data collected from the MRI machine
  %   ( kx, ky, slice, coil )
  % cf - the cutoff spatial frequency
  %
  % Optional Inputs:
  % invN - the inverse covariance matrix specifying how the noise between coils is correlated
  %
  % Outputs:
  % senseMaps - the sensitivity maps
  %
  % Written by Nicholas Dwork - Copyright 2019
  %
  % https://github.com/ndwork/dworkLib.git
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  p = inputParser;
  p.addRequired( 'kData', @isnumeric );
  p.addRequired( 'kcf', @ispositive );
  p.addParameter( 'doCheckAdjoint', false, @islogical );
  p.addParameter( 'initialGuess', [], @isnumeric );
  p.addParameter( 'invN', [], @isnumeric );
  p.addParameter( 'maxIterLSQR', 20, @ispositive );
  p.parse( kData, kcf, varargin{:} );
  doCheckAdjoint = p.Results.doCheckAdjoint;
  initialGuess = p.Results.initialGuess;
  invN = p.Results.invN;
  maxIterLSQR = p.Results.maxIterLSQR;

  [ Nx, Ny, nSlices, nCoils ] = size( kData );
  ks = size2fftCoordinates( [ Nx Ny ] );
  [ kxs, kys ] = meshgrid( ks{1}, ks{2} );
  kMag = sqrt( kxs .* kxs + kys .* kys );
  kMask = repmat( kMag < kcf, [ 1 1 nCoils ] );

  nSamples = Ny * Nx;

  function out = Z( x )
    out = zeros( Nx, Ny, nCoils );
    out( kMask == 1 ) = x(:);
  end

  function out = Zadj( x )
    out = x( kMask == 1 );
  end

  function out = F( x )
    tmp = ifftshift( ifftshift( x, 1 ), 2 );   % iShiftedX
    out = 1 / sqrt( nSamples ) * fft( fft( tmp, [], 1 ), [], 2 );
    out = fftshift( fftshift( out, 1 ), 2 );
  end

  function out = Fadj( y )
    tmp = ifftshift( ifftshift( y, 1 ), 2 );
    out = sqrt( nSamples ) * ifft( ifft( tmp, [], 1 ), [], 2 );
    out = fftshift( fftshift( out, 1 ), 2 );
  end

  function out = applyA( in, type )
    if nargin < 2, type = 'notransp'; end

    if strcmp( type, 'notransp' )
      Zs = Z( in );
      FadjZs = Fadj( Zs );
      XFadjZs = bsxfun( @times, FadjZs, recon );
      FXFadjZs = F( XFadjZs );
      out = FXFadjZs(:);

    else

      y = reshape( in, [ Nx, Ny, nCoils ] );
      Fadjy = Fadj( y );
      XadjFadjy = bsxfun( @times, Fadjy, conj( recon ) );
      FXadjFadjy = F( XadjFadjy );
      ZadjFXadjFadjy = Zadj( FXadjFadjy );
      out = ZadjFXadjFadjy(:);

    end
  end

  if doCheckAdjoint
    if ~checkAdjoint( kMask, @F, @Fadj ), error('checkAdjoint failed'); end
    if ~checkAdjoint( kMask(kMask==1), @Z, @Zadj ), error('checkAdjoint failed'); end
    if ~checkAdjoint( kMask(kMask==1), @applyA ), error('checkAdjoint failed'); end
  end

  fftSenseMaps0 = ones( sum( kMask(:) ), nSlices );
  if numel( initialGuess ) > 0
    fftInitialGuess = F( initialGuess );
    for sliceIndx = 1 : nSlices
      tmp = squeeze( fftInitialGuess( :, :, sliceIndx, : ) );
      fftSenseMaps0(:,sliceIndx) = tmp( kMask==1 );
    end
  end

  senseMaps = zeros( Nx, Ny, nSlices, nCoils );
  for sliceIndx = 1 : nSlices

    sliceKData = squeeze( kData( :, :, sliceIndx, : ) );

    fftSliceSenseMaps0 = fftSenseMaps0( :, sliceIndx );

    %Ax0 = applyA( fftSliceSenseMaps0(:), 'notransp' );
    %disp([ 'Error0: ', num2str( norm( Ax0(:) - sliceKData(:) ) ) ]);

    %TODO:  Provide warm starting with last best guess
    x = lsqr( @applyA, sliceKData(:), [], maxIterLSQR, [], [], fftSliceSenseMaps0 );
    sliceSenseMap = reshape( Fadj( Z( x ) ), [ Nx Ny nCoils ] );
    senseMaps(:,:,sliceIndx,:) = sliceSenseMap;

    %Ax = applyA( x(:), 'notransp' );
    %disp([ 'Error: ', num2str( norm( Ax(:) - sliceKData(:) ) ) ]);
    %showImageCube( sliceSenseMap, 'border', 2, 'borderValue', 'max' );

  end

end

