
function senseMaps = mccs_makeInitialSensitivityMap( kData, varargin )

  p = inputParser;
  p.addOptional( 'noiseCoords', [], @isnumeric );
  p.addOptional( 'kcf', [], @isnumeric );
  p.parse( varargin{:} );
  noiseCoords = p.Results.noiseCoords;
  kcf = p.Results.kcf;

  sKData = size( kData );
  nSamples = sKData(1) * sKData(2);

  shiftedRecons = sqrt( nSamples ) * ifftc( ifftc( kData, 1 ), 2 );
  coilRecons = fftshift( fftshift( shiftedRecons, 1 ), 2 );
  nSlices = size( coilRecons, 3 );
  nCoils = size( coilRecons, 4 );

  ssqRecon = norms( coilRecons, 2, 4 );

  if numel( kcf ) > 0

    ks = size2fftCoordinates( size(ssqRecon) );
    [ kxs, kys ] = meshgrid( ks{2}, ks{1} );
    kMag = sqrt( kxs .* kxs + kys .* kys );
    kcMask = kMag < max(kcf(:)) * 10;
    filteredCoilImgs = zeros( size( coilRecons ) );
    for sliceIndx = 1 : nSlices
      for coilIndx = 1 : nCoils
        fftCoilImg = fftshift( fft2( ifftshift( coilRecons(:,:,sliceIndx,coilIndx ) ) ) );
        fftFilteredCoilImgs = fftCoilImg .* kcMask;
        filteredCoilImgs(:,:,sliceIndx,coilIndx) = ...
          fftshift( ifft2( ifftshift( fftFilteredCoilImgs ) ) );
      end
      fftSsqCoilRecon = fftshift( fft2( ifftshift( ssqRecon ) ) );
      fftFilteredSsqCoilRecon = fftSsqCoilRecon .* kcMask;
      filteredCoilSsqRecon = fftshift( ifft2( ifftshift( fftFilteredSsqCoilRecon ) ) );
    end

  else

    filteredCoilImgs = coilRecons;
    filteredCoilSsqRecon = ssqRecon;

  end

  senseMaps = bsxfun( @rdivide, filteredCoilImgs, filteredCoilSsqRecon );

  %sMaps = size( senseMaps );
  %senseMaps = padData( senseMaps, [ 2*sMaps(1:2) size(senseMaps,3) size(senseMaps,4) ], ...
  %  mean( senseMaps(:) ) );


%  sMaps = size( senseMaps );
%   ks = size2fftCoordinates( sMaps(1:2) );
%   [ kxs, kys ] = meshgrid( ks{2}, ks{1} );
%   kMag = sqrt( kxs .* kxs + kys .* kys );
%   nCoils = size( senseMaps, 4 );
%   kcMask = kMag < kcf;
%   for sliceIndx = 1 : size( senseMaps, 3 )
%     for coilIndx = 1 : nCoils
%       fftMap = fftshift( fft2( ifftshift( senseMaps(:,:,sliceIndx,coilIndx ) ) ) );
%       filteredMap = fftMap .* kcMask;
%       senseMaps(:,:,sliceIndx,coilIndx) = fftshift( ifft2( ifftshift( filteredMap ) ) );
%     end
%   end


%   fftRecons = mri_fftRecon( kData );
%   noiseRegions = fftRecons( noiseCoords(2):noiseCoords(4), ...
%     noiseCoords(1):noiseCoords(3), :, : );
%   
%   noiseMean = mean( abs( noiseRegions(:) ) );
%   noiseStd = std( abs( noiseRegions(:) ) );
%   thresh = noiseMean + 3 * noiseStd;
%   signalMask = abs( fftRecons ) > thresh;
% 
%   senseMaps = noisyMaps .* signalMask;



  fitPlane = 1;
  if fitPlane == 1
    roughMaps = senseMaps;
    sMaps = size( roughMaps );
    imgCoords = size2imgCoordinates( 2 * sMaps(1:2) );
    [xs,ys] = meshgrid( imgCoords{2}, imgCoords{1} );
    cropxs = cropData( xs, size(ssqRecon) );
    cropys = cropData( ys, size(ssqRecon) );

    nCoils = size( roughMaps, 4 );
    absRecon = abs( ssqRecon );
    senseMaps = cell( 1, 1, 1, nCoils );
    cMask = [ [ 1 1 ]; [ 1 0 ]; ];
    %parfor coil = 1 : nCoils
    for coil = 1 : nCoils
      coilRoughMap = roughMaps(:,:,1,coil);
      realPolyCoefs = polyFit2( cropxs, cropys, real(coilRoughMap), 1, 1, 'w', absRecon, ...
        'cMask', cMask );
      imagPolyCoefs = polyFit2( cropxs, cropys, imag(coilRoughMap), 1, 1, 'w', absRecon, ...
        'cMask', cMask );
      senseMaps{coil} = evaluatePoly2( realPolyCoefs, xs, ys ) + ...
        1i * evaluatePoly2( imagPolyCoefs, xs, ys );
    end
    senseMaps = cell2mat( senseMaps );
  end

end

