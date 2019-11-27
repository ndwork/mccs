
function senseMaps = mccs_makeInitialSensitivityMap( kData )

  shiftedKData = ifftshift( ifftshift( kData, 1 ), 2 );
  
  sKData = size( kData );
  nSamples = sKData(1) * sKData(2);
  shiftedRecons = sqrt( nSamples ) * ifft( ifft( shiftedKData, [], 1 ), [], 2 );

  coilRecons = fftshift( fftshift( shiftedRecons, 1 ), 2 );

  ssqRecon = norms( coilRecons, 2, 4 );

  senseMaps = bsxfun( @rdivide, coilRecons, ssqRecon );

%   sMaps = size( roughMaps );
%   imgCoords = size2imgCoordinates( sMaps(1:2) );
%   [xs,ys] = meshgrid( imgCoords{1}, imgCoords{2} );
% 
%   nCoils = size( roughMaps, 4 );
%   absRecon = abs( ssqRecon );
%   senseMaps = cell( 1, 1, 1, nCoils );
%   %parfor coil = 1 : nCoils
% for coil = 1 : nCoils
%     coilRoughMap = roughMaps(:,:,1,coil);
%     polyCoefs = polyFit2( xs, ys, real(coilRoughMap), 1, 1, 'w', absRecon ) + ...
%                 1i * polyFit2( xs, ys, real(coilRoughMap), 1, 1, 'w', absRecon );
%     senseMaps{coil} = evaluatePoly2( polyCoefs, xs, ys );
%   end
%   senseMaps = cell2mat( senseMaps );
  
end

