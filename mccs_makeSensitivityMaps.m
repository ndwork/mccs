
function senseMaps = mccs_makeSensitivityMaps( recon, kData, kcf, lambda_s, varargin )
  % senseMaps = mccs_makeSensitivityMaps( recon, kData, kcf [, ...
  %   'maxIterOpt', maxIterOpt ] )
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
  p.addParameter( 'maxIterOpt', 100, @ispositive );
  p.addParameter( 'verbose', false, @islogical );
  p.parse( kData, kcf, varargin{:} );
  doCheckAdjoint = p.Results.doCheckAdjoint;
  initialGuess = p.Results.initialGuess;
  invN = p.Results.invN;
  maxIterOpt = p.Results.maxIterOpt;
  verbose = p.Results.verbose;

  [ Nx, Ny, nSlices, nCoils ] = size( kData );
  ks = size2fftCoordinates( [ Nx Ny ] );
  [ kxs, kys ] = meshgrid( ks{1}, ks{2} );
  kMag = sqrt( kxs .* kxs + kys .* kys );
  kBwMask = repmat( kMag < kcf, [ 1 1 nCoils ] );
  kMask = squeeze( sum( kData, 3 ) ) ~= 0;
  nPix = Ny * Nx;

  function out = Z( x )
    out = zeros( Nx, Ny, nCoils );
    out( kBwMask == 1 ) = x(:);
  end

  function out = Zadj( x )
    out = x( kBwMask == 1 );
  end

  function out = F( x )
    tmp = ifftshift( ifftshift( x, 1 ), 2 );   % iShiftedX
    out = 1 / sqrt( nPix ) * fft( fft( tmp, [], 1 ), [], 2 );
    out = fftshift( fftshift( out, 1 ), 2 );
  end

  function out = Fadj( y )
    tmp = ifftshift( ifftshift( y, 1 ), 2 );
    out = sqrt( nPix ) * ifft( ifft( tmp, [], 1 ), [], 2 );
    out = fftshift( fftshift( out, 1 ), 2 );
  end

  nb = sum( kMask(:) );
  function out = applyA( in, type )
    if nargin < 2, type = 'notransp'; end

    if strcmp( type, 'notransp' )
      Zs = Z( in );
      FadjZs = Fadj( Zs );
      XFadjZs = bsxfun( @times, FadjZs, recon );
      FXFadjZs = F( XFadjZs );
      DFXFadjZs = FXFadjZs( kMask == 1 );
      out = [ DFXFadjZs(:); FadjZs(:); FadjZs(:); ];

    else

      y1 = zeros( Nx, Ny, nCoils );
      y1( kMask == 1 ) = in( 1 : nb );
      Fadjy = Fadj( y1 );
      XadjFadjy = bsxfun( @times, Fadjy, conj( recon ) );
      FXadjFadjy = F( XadjFadjy );
      ZadjFXadjFadjy = Zadj( FXadjFadjy );
      out1 = ZadjFXadjFadjy(:);

      y2 = reshape( in( nb + 1 : nb + nPix*nCoils ), [ Nx, Ny, nCoils ] );
      Fy2 = F( y2 );
      ZadjFy2 = Zadj( Fy2 );
      out2 = ZadjFy2(:);

      y3 = reshape( in( nb + nPix*nCoils + 1 : end ), [ Nx, Ny, nCoils ] );
      Fy3 = F( y3 );
      ZadjFy3 = Zadj( Fy3 );
      out3 = ZadjFy3(:);

      out = out1 + out2 + out3;
    end
  end

  if doCheckAdjoint
    if ~checkAdjoint( kBwMask, @F, @Fadj ), error('checkAdjoint failed'); end
    if ~checkAdjoint( kBwMask(kBwMask==1), @Z, @Zadj ), error('checkAdjoint failed'); end
    if ~checkAdjoint( kBwMask(kBwMask==1), @applyA ), error('checkAdjoint failed'); end
  end

  fftSenseMaps0 = ones( sum( kBwMask(:) ), nSlices );
  if numel( initialGuess ) > 0
    fftInitialGuess = F( initialGuess );
    for sliceIndx = 1 : nSlices
      tmp = squeeze( fftInitialGuess( :, :, sliceIndx, : ) );
      fftSenseMaps0(:,sliceIndx) = tmp( kBwMask==1 );
    end
  end

  f = @(y) 0;
  proxf = @(x,t) x;

  function out = proxg1Conj( y, t )
    out = ( y - t * b ) ./ ( t + 1 );
  end

  function out = proxg2Conj( y, t )
    out = softThresh( y, t );
  end

  function out = proxg3Conj( y, lambda )
    S = reshape( y, [ nPix, nCoils ] );
    out = proxConjNucNorm( S, lambda );
    out = out( : );
  end

  proxgConj = @(y,t) [ proxg1Conj( y( 1 : nb ), t ); ...
                       proxg2Conj( y( nb+1 : nb+nPix*nCoils ), t ); ...
                       proxg3Conj( y( nb+nPix*nCoils+1 : end ), lambda_s ); ];

  function out = g( y )
    out = 0.5 * norm( y(1:nb) - b(:) ).^2;  % Note, doesn't include indicator function
    disp( 'NICK, YOU NEED TO ADD THIRD TERM HERE' );
  end

  senseMaps = zeros( Nx, Ny, nSlices, nCoils );
  for sliceIndx = 1 : nSlices

    sliceKData = squeeze( kData( :, :, sliceIndx, : ) );
    b = sliceKData( kMask == 1 );

    fftSliceSenseMaps0 = fftSenseMaps0( :, sliceIndx );

    %Ax0 = applyA( fftSliceSenseMaps0(:), 'notransp' );
    %disp([ 'Error0: ', num2str( norm( Ax0(:) - sliceKData(:) ) ) ]);

    beta = 1;
    tau = 1d2;  % Works well
    theta = 1.9;
    if verbose == true
      
      normA = powerIteration( @applyA, fftSliceSenseMaps0(:) );
      %[x,objValues] = chambollePock( fftSliceSenseMaps0(:), proxf, proxgConj, tau, ...
      %  'A', @applyA, 'normA', normA, 'N', maxIterOpt, 'verbose', verbose, ...
      %  'theta', theta, 'f', f, 'g', @g );   %#ok<ASGLU>
      [x,objValues] = chambollePockWLS( fftSliceSenseMaps0(:), proxf, proxgConj, ...
        'A', @applyA, 'beta', beta, 'N', maxIterOpt, 'verbose', verbose, ...
        'doCheckAdjoint', doCheckAdjoint, 'f', f, 'g', @g );   %#ok<ASGLU>

    else

      x = chambollePock( fftSliceSenseMaps0(:), proxf, proxgConj, tau, ...
        'A', @applyA, 'normA', normA, 'N', maxIterOpt, 'theta', theta );
      %x = chambollePockWLS( fftSliceSenseMaps0(:), proxf, proxgConj, 'A', @applyA, ...
      %  'beta', beta, 'N', maxIterOpt, 'verbose', verbose, 'doCheckAdjoint', doCheckAdjoint );

    end

    sliceSenseMap = reshape( Fadj( Z( x ) ), [ Nx Ny nCoils ] );
    senseMaps(:,:,sliceIndx,:) = sliceSenseMap;

    %Ax = applyA( x(:), 'notransp' );
    %disp([ 'Error: ', num2str( norm( Ax(:) - sliceKData(:) ) ) ]);
    %showImageCube( sliceSenseMap, 'border', 2, 'borderValue', 'max' );

  end

end


