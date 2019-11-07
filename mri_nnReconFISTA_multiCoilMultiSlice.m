
function recon = mri_nnReconFISTA_multiCoilMultiSlice( kData, senseMaps, lambda, varargin )
  % recon = mri_nnReconFISTA_multiCoilMultiSlice( kData, senseMaps, lambda [ , ...
  %   'checkAdjoint', true/false ] )
  %
  % Inputs:
  % kData - the k-space data collected from the MRI machine
  %   ( Nx, Ny, slice, coil )
  % senseMaps - an array of size ( Nx, Ny, nSlices, nCoils ) describing the 
  %   sensitivity of each coil
  % lambda - the FISTA regularization parameter
  %
  % Optional Inputs:
  % checkAdjoint - true/false
  % nIter - number of FISTA iterations
  %
  % Output:
  % recon - the reconstruction of size ( Nx, Ny, nSlices )
  %
  % Written by Nicholas Dwork - Copyright 2019
  %
  % https://github.com/ndwork/dworkLib.git
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular purpose.

  p = inputParser;
  p.addRequired( 'kData', @isnumeric );
  p.addRequired( 'senseMaps', @isnumeric );
  p.addRequired( 'lambda', @(x) x >= 0 );
  p.addParameter( 'checkAdjoints', false, @islogical );
  p.addParameter( 'debug', false, @(x) isnumeric(x) || islogical(x) );
  p.addParameter( 'initialGuess', [], @isnumeric );
  p.addParameter( 'nIter', [], @ispositive );
  p.addParameter( 'printEvery', 1, @ispositive );
  p.addParameter( 'verbose', false, @(x) isnumeric(x) || islogical(x) );
  p.parse( kData, senseMaps, lambda, varargin{:} );
  checkAdjoints = p.Results.checkAdjoints;
  debug = p.Results.debug;
  initialGuess = p.Results.initialGuess;
  nIter = p.Results.nIter;
  printEvery = p.Results.printEvery;
  verbose = p.Results.verbose;
  

  [ Nx, Ny, nSlices, ~ ] = size( kData );

  slices = cell( nSlices, 1 );
  parforObj = parforProgress( nSlices );
  parfor sliceIndx = 1 : nSlices
    parforObj.progress( sliceIndx );   %#ok<PFBNS>
    % TODO: parallelize this loop
    kData4Slice = squeeze( kData(:,:,sliceIndx,:) );
    senseMapOfSlice = squeeze( senseMaps(:,:,sliceIndx,:) );

    if numel( initialGuess ) > 0
      reconGuess = initialGuess(:,:,sliceIndx);
    else
      reconGuess = [];
    end

    thisRecon = nnReconFISTA_slice( kData4Slice, senseMapOfSlice, lambda, ...
      'reconGuess', reconGuess, 'nIter', nIter, 'printEvery', printEvery, ...
      'debug', debug, 'checkAdjoints', checkAdjoints, 'verbose', verbose );
    slices{sliceIndx} = thisRecon;
  end
  parforObj.clean;

  recon = zeros( Nx, Ny, nSlices );
  for sliceIndx = 1 : nSlices
    recon(:,:,sliceIndx) = slices{sliceIndx};
  end
end


function recon = nnReconFISTA_slice( samples, senseMaps, lambda, varargin )
  % recon = nnReconFISTA_slice( samples, senseMaps, lambda [, 'debug', debug, 'nIter', nIter, ...
  %   'polish', polish, 'printEvery', printEvery, 'verbose', verbose ] )
  %
  % This routine minimizes 0.5 * || A x - b ||_2^2 + lambda || X ||_*
  %   where A is sampleMask * Fourier Transform * SenseMaps * [ real part; imag part; ]..
  %
  % Inputs:
  % samples - a 2D array that is zero wherever a sample wasn't acquired
  % senseMaps - a 3D array specifying the sensitivity of each location for each coil
  % lambda - regularization parameter
  %
  % Optional Inputs:
  % debug - if true, reduces the default number of iterations to 30 and forces verbose
  %         statements during optimization
  % nIter - the number of iterations that FISTA will perform (default is 100)
  % polish - if set to true, adds a polishing step (default is false)
  % printEvery - FISTA prints a verbose statement every printEvery iterations
  % verbose - if true, prints informative statements
  %
  % Written by Nicholas Dwork - Copyright 2019
  %
  % https://github.com/ndwork/dworkLib.git
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular purpose.

  p = inputParser;
  p.addRequired( 'samples', @isnumeric );
  p.addRequired( 'senseMaps', @isnumeric );
  p.addRequired( 'lambda', @(x) x >= 0 );
  p.addParameter( 'checkAdjoints', false, @islogical );
  p.addParameter( 'debug', false, @(x) isnumeric(x) || islogical(x) );
  p.addParameter( 'nIter', [], @(x) ispositive(x) || numel(x) == 0 );
  p.addParameter( 'printEvery', 1, @ispositive );
  p.addParameter( 'reconGuess', [], @isnumeric );
  p.addParameter( 'verbose', false, @(x) isnumeric(x) || islogical(x) );
  p.parse( samples, senseMaps, lambda, varargin{:} );
  checkAdjoints = p.Results.checkAdjoints;
  debug = p.Results.debug;
  nIter = p.Results.nIter;
  printEvery = p.Results.printEvery;
  reconGuess = p.Results.reconGuess;
  verbose = p.Results.verbose;

  if numel( nIter ) == 0
    if debug == true
      nIter = 30;
    else
      nIter = 100;
    end
  end

  [ Nx, Ny, nCoils ] = size( samples );
  M = ( sum( abs(samples), 3 ) ~= 0 );
  nSamples = Nx * Ny;  % Note that A is square

  % S is the tall block matrix, where each block is the diagonal matrix of
  %   coil sensitivity values
  % A = M F S,  A' = S' F' M
  % A' * A = S' F' M F S
  % g = (1/2) || Ax - b ||_2^2
  % gGrad = A' A x - A' b;

  function out = F( x )
    Fx = 1/sqrt(nSamples) .* fftshift( fft2( ifftshift( ...
      x(:,:,1) + 1i * x(:,:,2) ) ) );
    out = cat( 3, real(Fx), imag(Fx) );
  end

  function out = Fadj( y )
    Fadjy = sqrt(nSamples) * ...
      fftshift( ifft2( ifftshift( y(:,:,1) + 1i * y(:,:,2) ) ) );
    out = cat( 3, real(Fadjy), imag(Fadjy) );
  end

  M2 = cat( 3, M, M );

  function out = A( x )
    out = zeros([ Nx Ny nCoils 2]);
    xComplex = x(:,:,1) + 1i * x(:,:,2);
    for coilIndx = 1 : nCoils
      Sx = senseMaps(:,:,coilIndx) .* xComplex;
      Sx = cat( 3, real(Sx), imag(Sx) );
      MFSx = M2 .* F( Sx );
      out(:,:,coilIndx,1) = MFSx(:,:,1);
      out(:,:,coilIndx,2) = MFSx(:,:,2);
    end
  end

  function out = Aadj( y )
    out = zeros([ Nx Ny ]);
    for coilIndx = 1 : nCoils
      My = M2 .* squeeze( y(:,:,coilIndx,:) );
      FadjMy = Fadj( My );
      FadjMy = FadjMy(:,:,1) + 1i * FadjMy(:,:,2);
      sAdjFadjMy = conj(senseMaps(:,:,coilIndx)) .* FadjMy;
      out = out + sAdjFadjMy;
    end
    out = cat( 3, real(out), imag(out) );
  end

  function out = applyA( x, type )
    if nargin < 2, type = 'notransp'; end

    if strcmp( type, 'notransp' )
      x = reshape( x, [Nx Ny 2] );
      out = A( x );
    else
      x = reshape( x, [Nx, Ny, nCoils, 2] );
      out = Aadj( x );
    end
    out = out(:);
  end

  b = cat( 4, real(samples), imag(samples) );
  function out = g( x )
    Ax = A( x );
    diff = Ax - b;
    out = 0.5 * norm( diff(:), 2 ).^2;
  end

  Aadjb = Aadj( b );
  function out = gGrad( x )
    out = Aadj( A( x ) ) - Aadjb;
  end

  if checkAdjoints == true
    % Variable used during debugging of this routine
    checkResult = nnReconFISTA_checkAdjoints( samples, @F, @Fadj, @A, @Aadj );
    if checkResult ~= false, disp([ 'Adjoints test passed' ]); end
  end

  function out = proxth( x, t )
    xImg = reshape( x(:,:,1), [ Nx Ny ] ) + 1i * reshape( x(:,:,2), [ Nx Ny ] );
    nn = proxNucNorm( xImg, t * lambda );
    out = cat( 3, real( nn ), imag( nn ) );
  end

  function out = h( x )
    xImg = reshape( x(:,:,1), [ Nx Ny ] ) + 1i * reshape( x(:,:,2), [ Nx Ny ] );
    out = lambda * nucNorm( xImg );
  end

  if numel( reconGuess ) == 0
    x0 = zeros( Nx, Ny, 2 );
    ssqRecon = mri_ssqRecon( samples );
    x0(:,:,1) = sqrt( ssqRecon );
  else
    x0 = cat( 3, real(reconGuess), imag(reconGuess) );
  end

  if lambda > 0
    t = 1;
    if debug
      %[recon,oValues] = fista( x0, @g, @gGrad, proxth, 'h', @h, 'verbose', verbose );   %#ok<ASGLU>
      [recon,oValues] = fista_wLS( x0, @g, @gGrad, proxth, 'h', @h, ...
        't0', t, 'N', nIter, 'verbose', true, 'printEvery', printEvery );   %#ok<ASGLU>
    else
      %recon = fista( x0, @g, @gGrad, proxth );   %#ok<UNRCH>
      recon = fista_wLS( x0, @g, @gGrad, @proxth, 't0', t, 'N', nIter, ...
        'verbose', verbose, 'printEvery', printEvery );
    end

  else

    recon = lsqr( @applyA, b(:), [], [], [], [], x0(:) );
    recon = reshape( recon, [Ny Nx 2] );

  end

  recon = recon(:,:,1) + 1i * recon(:,:,2);
end


function out = nnReconFISTA_checkAdjoints( samples, F, Fadj, A, Aadj )

  sSamples = size( samples );
  imgRand = rand( [ sSamples(1:2) 2 ] );

  if checkAdjoint( imgRand, F, Fadj ) ~= true
    error( 'Fadj is not the adjoint of F' );
  end

  if checkAdjoint( imgRand, A, Aadj ) ~= true
    error( 'Aadj is not the adjoint of A' );
  end

  out = true;
end

