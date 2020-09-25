
function bInvNoiseCov = make_bInvNoiseCov( nKData, noiseCov )
  nCoils = size( noiseCov, 1 );
  sCov = nKData * nCoils;

  invNoiseCov = inv( noiseCov );
  [~,s,~] = svd( invNoiseCov, 'econ' );
  invNoiseCov = invNoiseCov / s(1);

  rowIndxs = zeros( sCov, 1 );
  colIndxs = zeros( sCov, 1 );
  values = zeros( sCov, 1 );

  rowIndx = 1;
  dataIndx = 1;
  for bIndx = 1 : nKData
    for coilIndx = 1 : nCoils
      indxs = dataIndx : dataIndx + nCoils - 1;

      rowIndxs( indxs ) = rowIndx;
      colIndxs( indxs ) = (bIndx-1) * nCoils + 1 : bIndx * nCoils;
      values( indxs ) = invNoiseCov( coilIndx, : );

      rowIndx = rowIndx + 1;
      dataIndx = dataIndx + nCoils;
    end
  end

  bInvNoiseCov = sparse( rowIndxs, colIndxs, values );
    % At this point, bInvNoiseCov is the matrix where b is a vector where the data
    % for the second k-space point is concatenated onto the data for the first k-space, 
    % point, and the data for the third point is concatenated onto the result, etc.
    % But b is ordered so that the vector is (b^(1), b^(2), ..., b^(C)) where b^(i) is the
    % data of the i^th coil.  So we must permute the covariance matrix.

  origColOrder = reshape( 1 : sCov, [ nKData nCoils ] );
  newColOrder = origColOrder';
  newColOrder = newColOrder(:);
  I = speye( sCov );
  P = I( :, newColOrder );  % permutation matrix

  bInvNoiseCov = P * bInvNoiseCov * P';
  
  bInvNoiseCov = bInvNoiseCov ./ max( bInvNoiseCov(:) );
end

