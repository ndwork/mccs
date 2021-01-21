
function findSenseMapsBW

  sImg = [ 8192 8192 ];  % size of resulting maps
  res = 1.5d-3;  % resolution in meter per pixel

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
  xBWs = zeros( nCoils, 1 );
  yBWs = zeros( nCoils, 1 );
  for coilIndx = 1 : nCoils
    cw = coilWidths( coilIndx );
    cl = coilLengths( coilIndx );
    coil = [ [ -0.5 * cw; 0; -0.5 * cl; ] ...
              [  0.5 * cw; 0; -0.5 * cl; ] ...
              [  0.5 * cw; 0;  0.5 * cl; ] ...
              [ -0.5 * cw; 0;  0.5 * cl; ] ...
              [ -0.5 * cw; 0; -0.5 * cl; ] ];
    coil = coil';

    coords = size2fftCoordinates( sImg );
    [ xLocs, yLocs ] = meshgrid( sImg(2)*res * coords{2}, sImg(1)*res * coords{1} );
    locs = [ xLocs(:)  yLocs(:)  zeros(numel(yLocs),1) ];
    tmp = mri_computeSensitivityBiotSavart( coil, locs );
    bsSensitivity = reshape( tmp, sImg );
    bsSpectrum = fftshift( ufft2( bsSensitivity ) );
    fwtmMask = abs( bsSpectrum ) >= 0.1 * max( abs( bsSpectrum(:) ) );
    ks = size2fftCoordinates( size( fwtmMask ) );
    kys = ks{1};  kxs = ks{2};
    [ kxLocs, kyLocs ] = meshgrid( kxs, kys );
    yBW = max( kyLocs( fwtmMask > 0 ) ) / max( yLocs(:) );  % cycles / meter
    xBW = max( kxLocs( fwtmMask > 0 ) ) / max( xLocs(:) );  % cycles / meter
    
    xBWs(coilIndx) = xBW;
    yBWs(coilIndx) = yBW;
  end

  disp( xBWs );
  disp( yBWs );
end

