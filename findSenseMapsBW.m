
function findSenseMapsBW

  sImg = [ 2048 2048 ];  % size of resulting maps
  res = 1.5d-3;  % resolution in meter per pixel
  maskThresh = 0.1;

  % first coil is non-square from a birdcage coil according to 
  % "A fast and accurate simulator for the design of birdcage coils in MRI"
  % by Giovanetti et al.
  cl = 0.12;  % coil length in meters
  radius = 0.08;  % meters
  nBirdcageCoils = 8;
  cw = 2 * pi * radius / nBirdcageCoils;  % coil width in meters

  coilLengths(1) = cl;
  coilWidths(1) = cw;


	% second coil is a square surface coil with sizes according to
  % "32-Element Receiver-Coil Array for Cardiac Imaging" by Hardy et al.
  coilLengths(2) = 0.07;  % coil length in meters
  coilWidths(2) = 0.07;  % coil width in meters

  nCoils = numel( coilWidths );
  BWs = zeros( nCoils, 1 );  % bandwidths

  coords = size2imgCoordinates( sImg );
  [ xLocs, yLocs ] = meshgrid( res * coords{2}, res * coords{1} );
  locs = [ xLocs(:)  yLocs(:)  zeros(numel(yLocs),1) ];
  fov = res * sImg;
  kIndxs = size2imgCoordinates( sImg );
  kys = kIndxs{1} / fov(1);  kxs = kIndxs{2} / fov(2);
  [ kxLocs, kyLocs ] = meshgrid( kxs, kys );
  freqMags = sqrt( kxLocs.^2 + kyLocs.^2 );

  for coilIndx = 1 : nCoils
    cw = coilWidths( coilIndx );
    cl = coilLengths( coilIndx );
    coil = [ [ -0.5 * cw; 0; -0.5 * cl; ] ...
             [  0.5 * cw; 0; -0.5 * cl; ] ...
             [  0.5 * cw; 0;  0.5 * cl; ] ...
             [ -0.5 * cw; 0;  0.5 * cl; ] ...
             [ -0.5 * cw; 0; -0.5 * cl; ] ];
    coil = coil';

    tmp = mri_computeSensitivityBiotSavart( coil, locs );
    bsSensitivity = reshape( tmp, sImg );
    bsSpectrum = fftshift( ufft2( bsSensitivity ) );
    mtf = abs( bsSpectrum );
    freqMask = mtf >= maskThresh * max( mtf(:) );

    bw = max( freqMags( freqMask > 0 ) );
    BWs( coilIndx ) = bw;
  end

  disp( BWs );
end

