
function wavSplit = makeWavSplit( sImg )

  minSplit = 16;

  yTmp = sImg(1);
  ySplitSize = 1;
  while ( mod( yTmp, 1 ) == 0 )
    ySplitSize = ySplitSize * 2;
    yTmp = yTmp / 2;
  end
  ySplitSize = ySplitSize / 4;
  ySplitSize = min( ySplitSize, minSplit );

  xTmp = sImg(2);
  xSplitSize = 1;
  while ( mod( xTmp, 1 ) == 0 )
    xSplitSize = xSplitSize * 2;
    xTmp = xTmp / 2;
  end
  xSplitSize = xSplitSize / 4;
  xSplitSize = min( xSplitSize, minSplit );

  wavSplit = zeros( [ ySplitSize xSplitSize ] );
  wavSplit(1) = 1;
end

