
function run_mri_mccsRecon
  close all; clear; rng(2);

  datacase = 1;
  simCase = 1;
  doSimulation = false;
  checkAdjoints = false;
  lambda = 1d-18;
  verbose = true;

  if doSimulation == true
    [ kData, senseMaps ] = loadSimData( simCase );
  else
    %kData = loadDatacase( datacase );
kData = [];  load kData.mat;
kData = kData( :, :, 10, : );
  end

  

  %mask = mri_makeIntensityMask( kData );
  %if ~exist( 'senseMaps', 'var' )
  %  senseMaps = mri_makeSensitivityMaps( kData, 'mask', mask, 'verbose', true );
  %save( 'senseMaps.mat', 'senseMaps' );
  %load 'senseMaps.mat';
  %end

  %ssqRecon = mri_ssqRecon( kData );
  %profile on;
  %recon = mri_csReconFISTA_multiCoilMultiSlice( kData, senseMaps, lambda, ...
  %  'checkAdjoints', checkAdjoints, 'verbose', verbose );
  %profile off;
  %profile viewer;
  %figure; imshowscale( ssqRecon, 5 );  title( 'ssqRecon' );
  %figure; imshowscale( abs( recon ), 5 );  title( 'abs(recon)' );


  [recon,senseMaps] = mri_mccsRecon( kData, lambda, 'verbose', verbose );

  figure; imshowscale( recon, 5 );
  figure; showImageCube( senseMaps, 'border', 1, 'borderValue', 'max' );

end



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

