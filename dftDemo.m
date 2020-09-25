
function dftDemo
  close all; clear; rng(1);

  N = 128;
  n1 = 8;
  n2 = 9;

  x1 = cos( 2*pi * (0:N-1)' * n1 / N );
  fftPlotThis( x1 )

  x2 = cos( 2*pi * (0:N-1)' * n2 / N );
  fftPlotThis( x2 )

  xInBetween = cos( 2*pi * (0:N-1)' * 0.5*(n1+n2) / N );
  fftPlotThis( xInBetween );

end


function fftPlotThis( x )
  fftx1 = fftshift( ufft( x ) );

  figure;
  subplot( 2, 1, 1 );  plotnice( x );
  subplot( 2, 1, 2 );  stemnice( abs(fftx1) );
  ylim([ 0 80 ]);
end
