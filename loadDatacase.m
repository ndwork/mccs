
function data = loadDatacase( datacase )
  % Returns data hybrid state

  dataDir = '/Users/nicholasdwork/.DellEMS/Vols/disk2s1';

  switch datacase

    case 1
      data = readOldMriDataOrgData( [ dataDir, '/mriData.Org/P14/kspace' ] );
      data = ifft( data, [], 3 );
      data = data( :, :, 5:10:end, : );
      data = data( :, :, 4, : );

    case 2
      data = readOldMriDataOrgData( [ dataDir, '/mriData.Org/P17/kspace' ] );
      data = ifft( data, [], 3 );
      data = data( :, :, 5:10:end, : );
      data = data( :, :, 4, : );

    case 3
      [data,header] = read_MR_rawdata( [ dataDir, '/mccs/P54784.7' ] );   %#ok<ASGLU>
      data = squeeze( data );
      data = ifft( ifftshift( data, 3 ), [], 3 );
      data = squeeze( data(:,:,62,:) );

    otherwise
      error( 'This datacase doesn''t exist' );
  end

end
