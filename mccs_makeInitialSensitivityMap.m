
function senseMaps = mccs_makeInitialSensitivityMap( kData )

  shiftedKData = ifftshift( ifftshift( kData, 1 ), 2 );
  
  sKData = size( kData );
  nSamples = sKData(1) * sKData(2);
  shiftedRecons = sqrt( nSamples ) * ifft( ifft( shiftedKData, [], 1 ), [], 2 );

  coilRecons = fftshift( fftshift( shiftedRecons, 1 ), 2 );

  ssqRecon = norms( coilRecons, 2, 4 );

  senseMaps = bsxfun( @rdivide, coilRecons, ssqRecon );

end

