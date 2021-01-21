
function senseMaps = mri_iSENSEnn_makeSensitivityMaps( recon, kData, lambda_s, varargin )
  % senseMaps = mri_iSENSEnn_makeSensitivityMaps_bad( recon, kData [, ...
  %   'maxIterOpt', maxIterOpt ] )
  %
  % Use nuclear norm regularization on each sensitivity coil map
  %
  % Inputs:
  % recon - the current reconstructed image
  % kData - the k-space data collecconjSliceReconted from the MRI machine
  %   ( ky, kx, slice, coil )
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
  p.addParameter( 'debug', false, @islogical );
  p.addParameter( 'checkAdjoints', false, @islogical );
  p.addParameter( 'initialGuess', [], @isnumeric );
  p.addParameter( 'maxIterOpt', 100, @ispositive );
  p.addParameter( 'verbose', false, @islogical );
  p.parse( kData, varargin{:} );
  debug = p.Results.debug;
  checkAdjoints = p.Results.checkAdjoints;
  initialGuess = p.Results.initialGuess;
  maxIterOpt = p.Results.maxIterOpt;
  verbose = p.Results.verbose;

  [ Ny, Nx, nSlices, nCoils ] = size( kData );
  kMask = squeeze( sum( abs(kData), 3 ) ) ~= 0;
  nPix = Nx * Ny;

  if numel( initialGuess ) ~= 0
    senseMaps0 = permute( initialGuess, [1 2 4 3] );
  else
    senseMaps0 = ones( Ny, Nx, nCoils, nSlices );
  end

  if checkAdjoints == true
    tmp = rand( Ny, Nx, nCoils );
    sliceRecon = rand( Ny, Nx ) + 1i * rand( Ny, Nx );
    conjSliceRecon = conj( sliceRecon );
    if ~checkAdjoint( tmp, @F, @Fadj ), error( 'F and Fadj are not adjoints' ); end    
    if ~checkAdjoint( tmp, @X, @Xadj ), error( 'X and Xadj are not adjoints' ); end
    if ~checkAdjoint( tmp, @R, @Radj ), error( 'R and Radj are not adjoints' ); end
    if ~checkAdjoint( tmp, @A, @Aadj ), error( 'A and Aadj are not adjoints' ); end
  end

  senseMaps = zeros( Ny, Nx, nSlices, nCoils );
  t = 1;
  for sliceIndx = 1 : nSlices
    sliceKData = squeeze( kData( :, :, sliceIndx, : ) );
    b = sliceKData( kMask == 1 );

    sliceRecon = recon(:,:,sliceIndx);
    conjSliceRecon = conj( sliceRecon );
    Aadjb = Aadj( b );
    sliceSenseMaps0 = senseMaps0( :, :, :, sliceIndx );

    printEvery = 1;
    if debug == true
      %[x,oValues] = fista( sliceSenseMaps0, @g, @gGrad, proxth, 'h', @h, ...
      %  'N', maxIterOpt, 'verbose', verbose );   %#ok<ASGLU>
      [x,oValues] = fista_wLS( sliceSenseMaps0, @g, @gGrad, proxth, 'h', @h, ...
        't0', t, 'N', maxIterOpt, 'verbose', true, 'printEvery', printEvery );   %#ok<ASGLU>
    else
      %x = fista( x0, @g, @gGrad, proxth, 'N', maxIterOpt );   %#ok<UNRCH>
      x = fista_wLS( sliceSenseMaps0, @g, @gGrad, @proxth, 't0', t, 'N', maxIterOpt, ...
        'verbose', verbose, 'printEvery', printEvery );
    end

    sliceSenseMap = reshape( x, [ Ny Nx nCoils ] );
    senseMaps(:,:,sliceIndx,:) = sliceSenseMap;
  end

  function out = X( in )
    out = bsxfun( @times, in, sliceRecon );
  end

  function out = Xadj( in )
    out = bsxfun( @times, in, conjSliceRecon );
  end

  function out = F( in )
    iShiftedX = ifftshift( ifftshift( in, 1 ), 2 );
    out = 1 / sqrt( nPix ) * fft( fft( iShiftedX, [], 1 ), [], 2 );
    out = fftshift( fftshift( out, 1 ), 2 );
  end

  function out = Fadj( in )
    iShiftedIn = ifftshift( ifftshift( in, 1 ), 2 );
    out = sqrt( nPix ) * ifft( ifft( iShiftedIn, [], 1 ), [], 2 );
    out = fftshift( fftshift( out, 1 ), 2 );
  end

  function out = R( in )
    % Apply data mask
    out = in( kMask == 1 );
  end

  function out = Radj( in )
    out = zeros( Ny, Nx, nCoils );
    out( kMask == 1 ) = in(:);
  end

  function out = A( in )
    out = R( F( X( in ) ) );
  end

  function out = Aadj( in )
    out = Xadj( Fadj( Radj( in ) ) );
  end

  function out = g( s )
    S = reshape( s, [ Ny Nx nCoils ] );
    AS = A( S );
    out = 0.5 * norm( AS(:) - b, 2 )^2;
  end

  function out = gGrad( x )
    X = reshape( x, [ Ny Nx nCoils ] );
    out = Aadj( A( X ) ) - Aadjb;
  end

  function out = h( s )
    out = 0;
    for coil = 1 : nCoils
      out = out + nucNorm( s(:,:,coil) );
    end
    out = lambda_s * out;
  end

  function out = proxth( in, t )
    out = zeros( Ny, Nx, nCoils );
    for coil = 1 : nCoils
      out(:,:,coil) = proxNucNorm( in(:,:,coil), t * lambda_s );
    end
  end
end
