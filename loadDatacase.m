
function data = loadDatacase( datacase )

  dataDir = '/Users/nicholasdwork/.DellEMS/Vols/disk2s1/mriData.Org';

  switch datacase

    case 1
      data = readOldMriDataOrgData( [ dataDir, '/P14/kspace' ] );

    case 2
      data = readOldMriDataOrgData( [ dataDir, '/P17/kspace' ] );

    otherwise
      error( 'This datacase doesn''t exist' );
  end

end
