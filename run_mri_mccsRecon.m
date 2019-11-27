
function run_mri_mccsRecon
  close all; clear; rng(1);

  %datacases = [ 3 1 7 2 5 8 6 ];
datacases = 3;

  doCheckAdjoint = false;
  simCase = 1;
  doSimulation = false;
  %vdSigma = 0.25;
vdSigma = 0.3;
  maxOuterIter = 100;
  useNoiseCov = true;

  %sampleFractions = [ 0.1 0.15 0.2 0.25 0.3 0.4 ];
sampleFractions = 0.3;

  %lambda_x = 5d-4;  % 1d-5 works well without noise covariance
  %lambda_s = 1d-4;

  labmda_iSENSEnn_x = 1d-3;
  labmda_iSENSEnn_s = 1d-2;
  verbose = true;


  for datacase = datacases
    if doSimulation == true
      [ sub_kData, senseMaps ] = loadSimData( simCase );  %#ok<ASGLU>
    else
      dataMatFile = [ 'sub_kData_', indx2str(datacase,99), '.mat' ];
      if exist( dataMatFile, 'file' )
        load( dataMatFile, 'sub_kData', 'noiseCoords' );
      else
        [sub_kData,noiseCoords] = loadDatacase( datacase );
        save( dataMatFile, 'sub_kData', 'noiseCoords' );
      end
    end
    sub_kData = sub_kData ./ max( abs( sub_kData(:) ) );

    if datacase == 3
      %kcf = 0.008;
      kcf = 0.015;
      lambda_xs = 1d-10 * ones( numel( sampleFractions ), 1 );
      lambda_ss = 1d-9 * ones( numel( sampleFractions ), 1 );
      lambda_hs = 1d2 * ones( numel( sampleFractions ), 1 );
    else
      kcf = 0.005;
      lambda_xs = 1d-11 * ones( numel( sampleFractions ), 1 );
      lambda_ss = 1d-3 * ones( numel( sampleFractions ), 1 );
      lambda_hs = 1d2 * ones( numel( sampleFractions ), 1 );
    end

    datacaseDir = [ './out/datacase_', indx2str(datacase,100) ];
    if ~exist( datacaseDir, 'dir' ), mkdir( datacaseDir ); end
    logFile = [ datacaseDir, '/metrics.log' ];
    logID = fopen( logFile, 'w' );
    fprintf( logID, 'datacase, sampleFraction, nSamples, type, mdm, mdm-2, niqe, piqe\n' );
    fclose( logID );

    mask = mri_makeIntensityMask( sub_kData );
    fullySampledSenseMapsMatFile = ['fullySampledSenseMaps_', indx2str(datacase,99), '.mat' ];
    if exist( fullySampledSenseMapsMatFile, 'file' )
      load( fullySampledSenseMapsMatFile, 'fullySampledSenseMaps' );
    else
      fullySampledSenseMaps = mri_makeSensitivityMaps( sub_kData, 'mask', mask, 'verbose', verbose );
      save( fullySampledSenseMapsMatFile, 'fullySampledSenseMaps' );
    end

    fftRecon_fullySampled = mri_fftRecon( sub_kData, 'multiSlice', true );
    if useNoiseCov == true
      noiseRegions = fftRecon_fullySampled( noiseCoords(2):noiseCoords(4), ...
        noiseCoords(1):noiseCoords(3), :, : );
      sNoiseRegions = size( noiseRegions );
      noiseVecs = reshape( noiseRegions, [ prod(sNoiseRegions(1:3)) sNoiseRegions(4) ] );
      noiseCov = cov( noiseVecs );
      noiseCov = noiseCov ./ max( abs( noiseCov(:) ) );
    else
      noiseCov = [];
    end

    ssqRecon_fullySampled = mri_ssqRecon( sub_kData, 'multiSlice', true);
    printMetrics( logFile, datacase, 1.0, numel(mask), ssqRecon_fullySampled, ...
      'ssqFullySampled', datacaseDir, 'senseMaps', fullySampledSenseMaps );

    fullySampledSakeL1EspiritRecon = sakeRecon2D( sub_kData, 'type', 'espiritL1' );
    printMetrics( logFile, datacase, 1.0, numel(mask), fullySampledSakeL1EspiritRecon, ...
      'sakeL1EspiritRecon', datacaseDir, ssqRecon_fullySampled );

    nSamples = sum( mask(:) );
    %parfor sampleFractionIndx = 1 : numel( sampleFractions )
for sampleFractionIndx = 1 : numel( sampleFractions )
      sampleFraction = sampleFractions( sampleFractionIndx );
      lambda_x = lambda_xs( sampleFractionIndx );
      lambda_s = lambda_ss( sampleFractionIndx );
      lambda_h = lambda_hs( sampleFractionIndx );

      outDir = [ datacaseDir, '/', num2str(sampleFraction,'%3.2f') ];
      if ~exist( outDir, 'dir' ), mkdir( outDir ); end

      [ Nx, Ny, nSlices, nCoils ] = size( sub_kData );
      sMask = [ Nx Ny ];
      if sampleFraction == 1
        vdMask = ones( sMask );
      else
        sampleGoal = round( sampleFraction*prod(sMask) );
        vdMask = vdSampleMask( sMask, vdSigma*max(Nx,Ny), sampleGoal );
      end
      kMask = sum( sum( abs( sub_kData ), 3), 4 ) ~= 0;
      vdMask = vdMask & kMask;
      maskFig = figure;  imshowscale( vdMask, 5 );
      saveas( maskFig, [outDir, '/sampleMask.png'] );
      close( maskFig );

      kData = sub_kData;
      for slice=1:nSlices
        for coil=1:nCoils
          kData(:,:,slice,coil) = vdMask .* sub_kData(:,:,slice,coil);
        end
      end

      ssqRecon = mri_ssqRecon( kData );
      printMetrics( logFile, datacase, sampleFraction, nSamples, ssqRecon, 'undersampledSsq', ...
        outDir, ssqRecon_fullySampled );

      [mccsRecon,mccsMaps] = mri_mccsRecon( kData, lambda_x, lambda_s, lambda_h, 'kcf', kcf, ...
        'maxOuterIter', maxOuterIter, 'noiseCov', noiseCov, 'verbose', verbose, ...
        'doCheckAdjoint', doCheckAdjoint, 'outDir', [ outDir, '/mccs/' ] );
      printMetrics( logFile, datacase, sampleFraction, nSamples, mccsRecon, 'mccsRecon', ...
        outDir, ssqRecon_fullySampled, 'senseMaps', mccsMaps );

      [csSenseRecon,senseMaps] = mri_csSenseRecon( kData, labmda_iSENSEnn_x, 'verbose', verbose );
      printMetrics( logFile, datacase, sampleFraction, nSamples, csSenseRecon, 'csSense', ...
        outDir, ssqRecon_fullySampled, 'senseMaps', senseMaps );

      [nnSenseRecon,nnMaps] = mri_nnSenseRecon( kData, labmda_iSENSEnn_x, 'senseMaps', senseMaps, ...
        'verbose', verbose );
      printMetrics( logFile, datacase, sampleFraction, nSamples, nnSenseRecon, 'nnSense', ...
        outDir, ssqRecon_fullySampled, 'senseMaps', nnMaps );

      [iSENSEnn_Recon,iSENSEnn_Maps] = mri_iSENSEnn( kData, labmda_iSENSEnn_x, labmda_iSENSEnn_s, ...
        'verbose', verbose, 'doCheckAdjoint', doCheckAdjoint, 'outDir', [ outDir, '/iSENSEnn/' ] );   %#ok<ASGLU>
      printMetrics( logFile, datacase, sampleFraction, nSamples, iSENSEnn_Recon, 'iSENSEnn_Recon', ...
        outDir, ssqRecon_fullySampled, 'senseMaps', iSENSEnn_Maps );

      sakeEspiritRecon = sakeRecon2D( kData, 'type', 'espirit' );
      printMetrics( logFile, datacase, sampleFraction, nSamples, sakeEspiritRecon, 'sakeEspiritRecon', ...
        outDir, ssqRecon_fullySampled);

      sakeL1EspiritRecon = sakeRecon2D( kData, 'type', 'espiritL1' );
      printMetrics( logFile, datacase, sampleFraction, nSamples, sakeL1EspiritRecon, 'sakeL1EspiritRecon', ...
        outDir, ssqRecon_fullySampled);
    end
  end
end

