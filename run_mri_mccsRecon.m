
function run_mri_mccsRecon
  close all; clear; rng(1);

  %datacases = [ 1 7 3 2 5 8 6 ];
datacases = 3;

  doCheckAdjoint = false;
  simCase = 1;
  doSimulation = false;
  %vdSigma = 0.25;
vdSigma = 0.3;
  maxOuterIter = 100;
  useNoiseCov = true;

  sampleFractions = [ 0.1 0.15 0.2 0.25 0.3 0.4 ];

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
      load( fullySampledSenseMapsMatFile, 'fullySampledSenseMaps' );  %#ok<LOAD>
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
    ssqRecon_fullySampledFig = figure;  imshowscale( abs( ssqRecon_fullySampled ), 5 );
    saveas( ssqRecon_fullySampledFig, [datacaseDir, '/recon_ssqFullySampled.png'] );
    close( ssqRecon_fullySampledFig );
    printMetrics( logFile, datacase, 1.0, numel(mask), ssqRecon_fullySampled, 'ssq Fully Sampled' );

    fullySampledSenseRecons = bsxfun( @times, fullySampledSenseMaps, ssqRecon_fullySampled );
    fullySampledSenseReconFig = figure;  showImageCube( abs(fullySampledSenseRecons), 5 );
    titlenice( 'Sense Images - Fully Sampled' );
    saveas( fullySampledSenseReconFig, [datacaseDir, '/recon_fullySampledSenseImages.png'] );
    close( fullySampledSenseReconFig );

    fullySampledSakeL1EspiritRecon = sakeRecon2D( sub_kData, 'type', 'espiritL1' );
    sake_ReconFig_fullySampledFig = figure;  imshowscale( abs( fullySampledSakeL1EspiritRecon ), 5 );
    title('sake L1-ESPIRiT Recon');
    saveas( sake_ReconFig_fullySampledFig, [datacaseDir, '/recon_sakeEspiritL1_fullySampled.png'] );
    close( sake_ReconFig_fullySampledFig );
    printMetrics( logFile, datacase, 1.0, numel(mask), fullySampledSakeL1EspiritRecon, ...
      'sakeL1EspiritRecon', ssqRecon_fullySampled );

    nSamples = sum( mask(:) );
    parfor sampleFractionIndx = 1 : numel( sampleFractions )
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
        vdMask = vdSampleMask( sMask, vdSigma*max(Nx,Ny), ...
          round( sampleFraction*prod(sMask) ) );
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
      undersampledSsqFig = figure;
      imshowscale( abs( ssqRecon ), 5 );  titlenice( 'Undersampled ssqRecon' );
      saveas( undersampledSsqFig, [outDir, '/recon_ssq.png'] );
      close( undersampledSsqFig );
      printMetrics( logFile, datacase, sampleFraction, nSamples, ssqRecon, 'ssq', ...
        ssqRecon_fullySampled );

      [mccsRecon,mccsMaps] = mri_mccsRecon( kData, lambda_x, lambda_s, lambda_h, ...
        'kcf', kcf, 'maxOuterIter', maxOuterIter, 'noiseCov', noiseCov, ...
        'verbose', verbose, 'doCheckAdjoint', doCheckAdjoint, ...
        'outDir', [ outDir, '/mccs/' ] );
      mccsReconFig = figure; imshowscale( abs(mccsRecon), 5 );  title('mccs Recon');
      saveas( mccsReconFig, [outDir, '/recon_mccs.png'] );
      close( mccsReconFig );
      mapsFig = figure;  showImageCube( abs( mccsMaps ), 5 );
      saveas( mapsFig, [outDir, '/senseMaps_mccs.png'] );
      close( mapsFig );
      printMetrics( logFile, datacase, sampleFraction, nSamples, mccsRecon, 'mccsRecon', ...
        ssqRecon_fullySampled );
  
      [csSenseRecon,senseMaps] = mri_csSenseRecon( kData, 100 * lambda_x, ...
        'verbose', verbose );
      %if ~exist( senseMapFile, 'file' ), save( senseMapFile, 'senseMaps' ); end
      csSenseFig = figure; imshowscale( abs(csSenseRecon), 5 );  title('csSense Recon');
      saveas( csSenseFig, [outDir, '/recon_csSense.png'] );
      close( csSenseFig );
      printMetrics( logFile, datacase, sampleFraction, nSamples, csSenseRecon, 'csSense', ...
        ssqRecon_fullySampled );

      nnSenseRecon = mri_nnSenseRecon( kData, lambda_x*1000, 'senseMaps', senseMaps, ...
        'verbose', verbose );
      nnSenseFig = figure; imshowscale( abs(nnSenseRecon), 5 );  title('nnSense Recon');
      saveas( nnSenseFig, [outDir, '/recon_nnSense.png'] );
      close( nnSenseFig );
      printMetrics( logFile, datacase, sampleFraction, nSamples, nnSenseRecon, 'nnSense', ...
        ssqRecon_fullySampled );
  
      [iSENSEnn_Recon,iSENSEnn_Maps] = mri_iSENSEnn( kData, labmda_iSENSEnn_x, labmda_iSENSEnn_s, ...
        'verbose', verbose, 'doCheckAdjoint', doCheckAdjoint, 'outDir', [ outDir, '/iSENSEnn/' ] );   %#ok<ASGLU>
      iSENSEnn_ReconFig = figure; imshowscale( abs(iSENSEnn_Recon), 5 );  title('mccs Recon');
      saveas( iSENSEnn_ReconFig, [outDir, '/recon_iSENSEnn.png'] );
      close( iSENSEnn_ReconFig );
      printMetrics( logFile, datacase, sampleFraction, nSamples, iSENSEnn_Recon, 'iSENSEnn_Recon', ...
        ssqRecon_fullySampled );

      sakeEspiritRecon = sakeRecon2D( kData, 'type', 'espirit' );
      sake_ReconFig = figure;  imshowscale( abs( sakeEspiritRecon ), 5 );  title('sake ESPIRiT Recon');
      saveas( sake_ReconFig, [outDir, '/recon_sakeEspirit.png'] );
      close( sake_ReconFig );
      printMetrics( logFile, datacase, sampleFraction, nSamples, sakeEspiritRecon, 'sakeEspiritRecon', ...
        ssqRecon_fullySampled);

      sakeL1EspiritRecon = sakeRecon2D( kData, 'type', 'espiritL1' );
      sake_ReconFig = figure;  imshowscale( abs( sakeL1EspiritRecon ), 5 );  title('sake L1-ESPIRiT Recon');
      saveas( sake_ReconFig, [outDir, '/recon_sakeEspiritL1.png'] );
      close( sake_ReconFig );
      printMetrics( logFile, datacase, sampleFraction, nSamples, sakeL1EspiritRecon, 'sakeL1EspiritRecon', ...
        ssqRecon_fullySampled);
    end
  end
end

