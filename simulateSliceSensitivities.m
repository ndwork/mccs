
function sensitivities = simulateSliceSensitivities( img, varargin )

  p = inputParser;
  p.addParameter( 'H', 0.5, @ispositive );
  p.addParameter( 'type', 'axial', @(x) true );
  p.parse( varargin{:} );
  H = p.Results.H;
  type = p.Results.type;

  sImg = size( img );
  D = sImg(1) / sImg(2) * H;  % depth of slice of interest

  coils = makeBirdcageCoils( );

  nCoils = numel( coils );

  sensitivities = zeros( [ sImg nCoils ] );

  coords = size2fftCoordinates( sImg );

  if strcmp( type, 'axial' )
    [ xs, ys ] = meshgrid( 0.5*H * coords{2}, 0.5*D * coords{1} );
    locs = [ xs(:)  ys(:)  zeros(numel(ys),1) ];
  elseif strcmp( type, 'sagittal' )
    [ ys, zs ] = meshgrid( 0.5*H * coords{2}, 0.5*D * coords{1} );
    locs = [ zeros(numel(ys),1)  ys(:)  zs(:)  ];
  end


  for coilIndx = 1 : nCoils
    thisCoil = coils{ coilIndx }';

    tmp = mri_computeSensitivityBiotSavart( thisCoil, locs );
    sensitivities( :, :, coilIndx ) = reshape( tmp, sImg );
  end

end

