
function [data,trueRecon] = loadSimData( simCase )

  switch simCase
    case 1
      img = rgb2gray( double( imread( '../data/17868sagittal1_64.png' ) ) / 255. );
      img = cropData( img, [ 176 256 ] );
      img = img / max( img(:) );
      trueRecon = img;
      sImg = size( img );

      sensitivities = simulateSliceSensitivities( img );
      sensitivities = sensitivities / max( sensitivities(:) );

      coilImgs = zeros( [ sImg, nCoils ] );
      fftCoilImgs = zeros( [ sImg, nCoils ] );
      for coilIndx = 1 : nCoils
        coilImgs(:,:,coilIndx) = sensitivities(:,:,coilIndx) .* img;
        fftCoilImgs(:,:,coilIndx) = fftshift( ufft2( ifftshift( coilImgs(:,:,coilIndx) ) ) );
      end

      %noise = 0.05 * ( rand( size( fftCoilImgs ) ) + 1i * rand( size( fftCoilImgs ) ) );
      noise = 0;
      data = fftCoilImgs + noise;
      data = reshape( data, [ sImg 1 nCoils ] );

    case 2
      img = rgb2gray( double( imread( '../data/27698transverse2_64.png' ) ) / 255. );
      img = cropData( img, [ 256 256 ] );
      img = img / max( img(:) );
      trueRecon = img;
      sImg = size( img );

      sensitivities = simulateSliceSensitivities( img, 'type', 'axial' );
      sensitivities = sensitivities / max( sensitivities(:) );

      nCoils = size( sensitivities, 3 );
      sensitivities = reshape( sensitivities, [ prod(sImg) nCoils ] );
      [u,s,v] = svd( sensitivities, 'econ' );
      s(8,8)=0;  s(7,7)=0;  s(6,6)=0;
      sensitivities = reshape( u * s * v', [ sImg nCoils ] );

      coilImgs = zeros( [ sImg, nCoils ] );
      fftCoilImgs = zeros( [ sImg, nCoils ] );
      for coilIndx = 1 : nCoils
        coilImgs(:,:,coilIndx) = sensitivities(:,:,coilIndx) .* img;
        fftCoilImgs(:,:,coilIndx) = fftshift( ufft2( ifftshift( coilImgs(:,:,coilIndx) ) ) );
      end

      %noise = 0.05 * ( rand( size( fftCoilImgs ) ) + 1i * rand( size( fftCoilImgs ) ) );
      noise = 0;
      data = fftCoilImgs + noise;
      data = reshape( data, [ sImg 1 nCoils ] );

    case 3
      img = rgb2gray( double( imread( 'bac.png' ) ) / 255. );
      img = imresize( img, [ 256 256 ], 'bicubic' );
      img = img / max( img(:) );
      trueRecon = img;
      sImg = size( img );
      nCoils = 4;
      coils = zeros( [ 2*sImg nCoils ] );
      d1 = sImg(1) / 2;
      d2 = sImg(2) / 2;

      radImg = makeRadialImg( 2*sImg, round( [ d1+sImg(1)/2 d2+1 ] ) );
      coils(:,:,1) = radImg / max( radImg(:) );
      radImg = makeRadialImg( 2*sImg, round( [ d1+sImg(1)/2 d2+sImg(2)/4 ] ) );
      coils(:,:,2) = radImg / max( radImg(:) );
      radImg = makeRadialImg( 2*sImg, round( [ d1+sImg(1)/2 d2+sImg(2)*3/4 ] ) );
      coils(:,:,3) = radImg / max( radImg(:) );
      radImg = makeRadialImg( 2*sImg, round( [ d1+sImg(1)/2 d2+sImg(2) ] ) );
      coils(:,:,4) = radImg / max( radImg(:) );

      %radImg = makeRadialImg( 2*sImg, round( [ d1+1 d2+sImg(2) ] ) );
      %coils(:,:,1) = radImg / max( radImg(:) );
      %radImg = makeRadialImg( 2*sImg, round( [ d1+sImg(1)/2 d2+sImg(2) ] ) );
      %coils(:,:,2) = radImg / max( radImg(:) );
      %radImg = makeRadialImg( 2*sImg, round( [ d1+sImg(1) d2+sImg(2) ] ) );
      %coils(:,:,3) = radImg / max( radImg(:) );
      %radImg = makeRadialImg( 2*sImg, round( [ d1+sImg(1)/2 d2+1 ] ) );
      %coils(:,:,4) = radImg / max( radImg(:) );
      coils = 1 - coils;
      coils = coils.^5;

      ks = size2fftCoordinates( 2 * sImg(1:2) );
      [ kxs, kys ] = meshgrid( ks{2}, ks{1} );
      kMag = sqrt( kxs .* kxs + kys .* kys );
      kcMask = kMag < kcf;

      extCoilMaps = zeros( [ 2*sImg, nCoils ] );
      for coilIndx = 1 : nCoils
        thisCoilMap = coils(:,:,coilIndx);
        fftCoilMap = fftshift( ufft2( ifftshift( thisCoilMap ) ) );
        maskedFftCoilMap = kcMask .* fftCoilMap;
        thisExtCoilMaps = fftshift( uifft2( ifftshift( maskedFftCoilMap ) ) );
        extCoilMaps(:,:,coilIndx) = thisExtCoilMaps / max( thisExtCoilMaps(:) );
      end

      extCoilCols = reshape( extCoilMaps, [ prod(2*sImg) nCoils ] );
      [u,s,v] = svd( extCoilCols, 'econ' );
      s(4,4)=0;
      extCoilMaps = reshape( u * s * v', [ 2*sImg nCoils ] );
      coilMaps = cropData( extCoilMaps, [ sImg nCoils ] );

      coilImgs = zeros( [ sImg, nCoils ] );
      fftCoilImgs = zeros( [ sImg, nCoils ] );
      for coilIndx = 1 : nCoils
        coilImgs(:,:,coilIndx) = coilMaps(:,:,coilIndx) .* img;
        fftCoilImgs(:,:,coilIndx) = fftshift( ufft2( ifftshift( coilImgs(:,:,coilIndx) ) ) );
      end

      %noise = 0.05 * ( rand( size( fftCoilImgs ) ) + 1i * rand( size( fftCoilImgs ) ) );
      noise = 0;
      data = fftCoilImgs + noise;
      data = reshape( data, [ sImg 1 nCoils ] );

    otherwise
      error( 'Unrecognized simCase' );
  end
  
  
end
