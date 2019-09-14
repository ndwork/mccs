
function [kData,senseMaps] = loadSimData( simCase )

  switch simCase
    case 1
      img = double( imread( 'cameraman.tif' ) ) / 255.;
      sImg = size( img );
      xs = ones(sImg(1),1) * (1:sImg(2));
      sMap1 = xs / max(xs(:));
      sMap2 = fliplr( sMap1 );
      ys = (1:sImg(1))' * ones(1,sImg(2));
      sMap3 = ys / max(ys(:));
      sMap4 = flipud( sMap3 );
      senseMaps = zeros([ sImg(1:2), 1, 4 ]);
      senseMaps(:,:,1,1) = sMap1;
      senseMaps(:,:,1,2) = sMap2;
      %senseMaps(:,:,1,3) = sMap3;
      %senseMaps(:,:,1,4) = sMap4;

      %senseMaps = 1i * senseMaps;

    otherwise
      error( 'Unrecognized simCase' );
  end

  sSenseMaps = size( senseMaps );
  kData = zeros( sSenseMaps );
  nSamples = sSenseMaps(1) * sSenseMaps(2);
  for slice = 1 : sSenseMaps(3)
    for coil = 1 : sSenseMaps(4)
      kData(:,:,slice,coil) = 1/sqrt(nSamples) * ...
        fftshift( fft2( ifftshift( senseMaps(:,:,slice,coil) .* img ) ) );
    end
  end

end
