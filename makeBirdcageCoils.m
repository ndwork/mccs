

function coils = makeBirdcageCoils( varargin )
  % coils = makeBirdcageCoils( varargin )
  %
  % Optional Inputs:
  % nCoils - number of coils; default is 8
  % R - radius, default is 0.2
  % L - the length of the birdcage coil
  %
  % Written by Nicholas Dwork, Copyright 2020

  p = inputParser;
  p.addParameter( 'nCoils', 8, @ispositive );
  p.addParameter( 'radius', 0.3, @ispositive );
  p.addParameter( 'length', 0.05, @ispositive );
  p.parse( varargin{:} );
  nCoils = p.Results.nCoils;
  radius = p.Results.radius;
  length = p.Results.length;

  coilAngles = 0 : 2 * pi / nCoils : 2*pi;
  coilAngles = coilAngles(1:nCoils) + coilAngles(2)/2;

  cw = 2 * pi * radius / nCoils;  % coil width

  coil = [ [ 0; -0.5 * cw; -0.5 * length; ] ...
           [ 0;  0.5 * cw; -0.5 * length; ] ...
           [ 0;  0.5 * cw;  0.5 * length; ] ...
           [ 0; -0.5 * cw;  0.5 * length; ] ...
           [ 0; -0.5 * cw; -0.5 * length; ] ];

  translation = [ radius; 0; 0; ];  % Translation
  coil = bsxfun( @plus, coil, translation );

  coils = cell( nCoils, 1 );
  for coilIndx = 1 : nCoils
    R = [ [ cos( coilAngles( coilIndx ) )  -sin( coilAngles( coilIndx ) )  0 ]; ...
          [ sin( coilAngles( coilIndx ) )   cos( coilAngles( coilIndx ) )  0 ]; ...
          [ 0 0 1 ] ];

    coils{ coilIndx } = R * coil;
  end

end

