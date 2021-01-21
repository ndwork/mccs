
function senseMaps = mccs_makeSensitivityMaps( recon, kData, kcf, lambda_s, ...
  lambda_h, varargin )
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
  % noiseCov - the nCoils x nCoils noise covariance matrix
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
  p.addParameter( 'maxIterOpt', 100, @ispositive );
  p.addParameter( 'noiseCov', [], @isnumeric );
  p.addParameter( 'verbose', false, @islogical );
  p.parse( kData, kcf, varargin{:} );
  doCheckAdjoint = p.Results.doCheckAdjoint;
  initialGuess = p.Results.initialGuess;
  maxIterOpt = p.Results.maxIterOpt;
  noiseCov = p.Results.noiseCov;
  verbose = p.Results.verbose;

  [ Ny, Nx, nSlices, nCoils ] = size( kData );
  nPix = Ny * Nx;
  kMask = squeeze( sum( abs(kData), 3 ) ) ~= 0;
  nb = sum( kMask(:) );
  kMaskCoil = sum( abs( kMask ), 3 ) > 0;
  nbPerCoil = sum( kMaskCoil(:) );

  ks = size2fftCoordinates( 2 * [ Ny Nx ] );
  [ kxs, kys ] = meshgrid( ks{2}, ks{1} );
  if numel( kcf ) == 1
    kMag = sqrt( kxs .* kxs + kys .* kys );
    kcMask = kMag < kcf;
  else
    kcMask = ( abs(kys) < kcf(1) ) & ( abs(kxs) < kcf(2) );
  end
  kcMask = repmat( kcMask, [ 1 1 nCoils ] );
  nHighFreq = sum( kcMask(:) == 0 );

  if numel( noiseCov ) > 0
    invNoiseCov = inv( noiseCov );
    [~,s,~] = svd( invNoiseCov, 'econ' );
    invNoiseCov = invNoiseCov ./ s(1);
    L = chol( invNoiseCov, 'lower' );
  end
  function out = applyL( in, type )
    if numel( noiseCov ) == 0, out = in; return; end

    % Assumes last dimension is coil dimension
    if nargin < 2, type = 'notransp'; end

    sIn = size( in );
    reshaped = reshape( in, [ prod( sIn(1:end-1) ) sIn(end) ] );
    if strcmp( type, 'notransp' )
      out = transpose( L * transpose( reshaped ) );
    else
      out = transpose( L' * transpose( reshaped ) );
    end
    out = reshape( out, sIn );
  end

  function out = F( x )
    iShiftedX = ifftshift( ifftshift( x, 1 ), 2 );
    out = 1 / sqrt( size(x,1) * size(x,2) ) * fft( fft( iShiftedX, [], 1 ), [], 2 );
    out = fftshift( fftshift( out, 1 ), 2 );
  end

  function out = Fadj( y )
    tmp = ifftshift( ifftshift( y, 1 ), 2 );
    out = sqrt( size(y,1) * size(y,2) ) * ifft( ifft( tmp, [], 1 ), [], 2 );
    out = fftshift( fftshift( out, 1 ), 2 );
  end

  function out = applyA1( in, type )
    if strcmp( type, 'notransp' )
      cropped = cropData( in, [ Ny Nx nCoils ] );
      Xin = bsxfun( @times, cropped, recon );
      FXin = F( Xin );
      DFXin = FXin( kMask == 1 );
      if numel( noiseCov ) > 0
        LStarDFXin = applyL( reshape( DFXin, [ nbPerCoil nCoils ] ), 'transp' );
      else
        LStarDFXin = DFXin;
      end
      out = LStarDFXin(:);
    else
      if numel( noiseCov ) > 0
        Lin = applyL( reshape( in, [ nbPerCoil nCoils ] ), 'notransp' );
      else
        Lin = in;
      end
      DTLin = zeros( Ny, Nx, nCoils );
      DTLin( kMask == 1 ) = Lin(:);
      FadjDTLin = Fadj( DTLin );
      XconjFadjDTLin = bsxfun( @times, FadjDTLin, conj( recon ) );
      out = padData( XconjFadjDTLin, [ 2*Ny 2*Nx nCoils ] );
    end
  end

%winY = kaiser( Ny, 5 );
%winX = kaiser( Nx, 5 );
%coilWin = winY' * winX;

  function out = applyA2( in, type )
    if strcmp( type, 'notransp' )
      s = reshape( in, [ 2*Ny 2*Nx nCoils ] );
      %Hs = bsxfun( @times, s, coilWin );
      Fs = F( s );
      DcFs = Fs( kcMask == 0 );
      out = sqrt( lambda_h ) * DcFs;
    else
      tmp = zeros( 2*Ny, 2*Nx, nCoils );
      tmp( kcMask == 0 ) = in;
      DbwIn = sqrt( lambda_h ) * tmp;
      out = Fadj( DbwIn );
      %out = bsxfun( @times, FadjDbwIn, coilWin );
    end
  end

  function out = applyA( in, type )
    if nargin < 2, type = 'notransp'; end

    if strcmp( type, 'notransp' )
      Ain1 = applyA1( reshape( in, [ 2*Ny 2*Nx nCoils ] ), type );
      Ain2 = applyA2( reshape( in, [ 2*Ny 2*Nx nCoils ] ), type );
      out = [ Ain1(:); Ain2(:); in(:); ];

    else
      y1 = in( 1 : nb );
      y2 = in( nb + 1 : nb + nHighFreq );
      y3 = in( nb + nHighFreq + 1 : nb + nHighFreq + 4*nPix*nCoils );

      A1T = applyA1( y1, type );
      A2T = applyA2( y2, type );

      out = A1T(:) + A2T(:) + y3(:);
    end
  end

  if doCheckAdjoint
    if ~checkAdjoint( kcMask, @F, @Fadj ), error('checkAdjoint failed'); end
    if ~checkAdjoint( kcMask, @applyA1 ), error('checkAdjoint failed'); end
    if ~checkAdjoint( kcMask, @applyA2 ), error('checkAdjoint failed'); end
    if ~checkAdjoint( kcMask, @applyA ), error('checkAdjoint failed'); end
  end

  f = @(s) lambda_s * nucNorm( reshape( s, [ 4*nPix nCoils ] ) );

  function out = g( y )
    y1 = y( 1 : nb );
    y2 = y( nb + 1 : nb + nHighFreq );
    y3 = y( nb + nHighFreq + 1 : nb + nHighFreq + 4*nPix*nCoils );  %#ok<NASGU>

    out1 = 0.5 * norm( y1 - LStarb(:) )^2;
    out2 = 0.5 * norm( y2 )^2;
    % I will not include the indicator function
    % I may want to indicate how much the constraint is violated.

    out = out1 + out2;
  end

  function out = proxf( y, t )
    out = proxNucNorm( reshape( y, [ 4*nPix nCoils ] ), t * lambda_s );
    out = out(:);
  end

  function out = proxgConj( y, t )
    y1 = y( 1 : nb );
    y2 = y( nb + 1 : nb + nHighFreq );
    y3 = y( nb + nHighFreq + 1 : nb + nHighFreq + 4*nPix*nCoils );

    out1 = ( y1 - t * LStarb(:) ) ./ ( t + 1 );
    out2 = y2 ./ ( t + 1 );    
    out3 = softThresh( y3, t );

    out = [ out1(:); out2(:); out3(:); ];
  end


  if numel( initialGuess ) == 0
    initialGuess = ones( [ 2*Ny 2*Nx nSlices nCoils ] );
  end

  senseMaps = initialGuess;
  for sliceIndx = 1 : nSlices
    sliceKData = squeeze( kData( :, :, sliceIndx, : ) );
    sliceGuess = initialGuess( :, :, sliceIndx, : );

    bKData = reshape( sliceKData, [nPix nCoils] );
    b = bKData( kMaskCoil(:) == 1, : );
    LStarb = applyL( b, 'transp' );

    beta = 1;
    if verbose == true
      [x,objValues] = pdhgWLS( sliceGuess(:), @proxf, @proxgConj, ...
        'A', @applyA, 'beta', beta, 'N', maxIterOpt, 'verbose', verbose, ...
        'doCheckAdjoint', doCheckAdjoint, 'f', f, 'g', @g );   %#ok<ASGLU>
    else
      x = pdhgWLS( sliceGuess(:), @proxf, @proxgConj, 'A', @applyA, ...
        'beta', beta, 'N', maxIterOpt, 'verbose', verbose, ...
        'doCheckAdjoint', doCheckAdjoint );
    end

    sliceSenseMap = reshape( x, [ 2*Ny 2*Nx nCoils ] );
    %senseMaps(:,:,sliceIndx,:) = cropData( sliceSenseMap, [ Ny Nx nCoils ] );
    senseMaps(:,:,sliceIndx,:) = sliceSenseMap;

    %Ax = applyA( x(:), 'notransp' );
    %disp([ 'Error: ', num2str( norm( Ax(:) - sliceKData(:) ) ) ]);
    %showImageCube( sliceSenseMap, 'border', 2, 'borderValue', 'max' );

  end

%figure;  plotnice( objValues );
end


