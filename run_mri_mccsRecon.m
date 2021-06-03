
function run_mri_mccsRecon
  close all; clear; rng(1);

  %datacases = [ 1, 0, 2, 6, 11, 5, 9, 3, 7, 8, 10 ];
datacases = [ 1 ];
  % Note, datacase 0 means simulation

  doCheckAdjoint = false;
  simCase = 2;
  vdSigma = 0.3;
  maxOuterIter = 50;

  sampleFractions = [ 0.10 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5 ];
  verbose = true;

  if ~exist( 'sub_kData', 'dir' ), mkdir( 'sub_kData' ); end
  if ~exist( 'fullySampledSenseMaps', 'dir' ), mkdir( 'fullySampledSenseMaps' ); end

  for datacase = datacases
    dataMatFile = [ 'sub_kData/sub_kData_', indx2str(datacase,99), '.mat' ];
    if exist( dataMatFile, 'file' )
      loadData = false;
      [ ~, noiseCoords, kcf, res, lambda_xs, lambda_ss, lambda_hs, lambda_iSENSEnn_x, ...
        lambda_iSENSEnn_s, lambda_csSENSE, lambda_nnSENSE, trueRecon ] = loadDatacase( ...
        datacase, sampleFractions, loadData, 'simCase', simCase );
      load( dataMatFile, 'sub_kData' );
    else
      loadData = true;
      [ sub_kData, noiseCoords, kcf, res, lambda_xs, lambda_ss, lambda_hs, lambda_iSENSEnn_x, ...
        lambda_iSENSEnn_s, lambda_csSENSE, lambda_nnSENSE, trueRecon ] = loadDatacase( ...
        datacase, sampleFractions, loadData, 'simCase', simCase );
      save( dataMatFile, 'sub_kData' );
    end

    if datacase == 0
      useNoiseCov = false;
    else
      useNoiseCov = true;
    end

    sub_kData = sub_kData ./ max( abs( sub_kData(:) ) );
    if numel( trueRecon ) > 0
      trueRecon = trueRecon ./ max( abs( sub_kData(:) ) );
    end

    datacaseDir = [ './out/datacase_', indx2str(datacase,100) ];
    if ~exist( datacaseDir, 'dir' ), mkdir( datacaseDir ); end
    logFile = [ datacaseDir, '/metrics.log' ];
    logID = fopen( logFile, 'w' );
    fprintf( logID, [ 'datacase, sampleFraction, nSamples, type, mdm, mdm-2, niqe, ', ...
      'piqe, mse, mae, correlation, angleErr, mi \n' ] );
    fclose( logID );

    mask = mri_makeIntensityMask( sub_kData );
    fullySampledSenseMapsMatFile = ['fullySampledSenseMaps/fullySampledSenseMaps_', ...
      indx2str(datacase,99), '.mat' ];
    if ~exist( fullySampledSenseMapsMatFile, 'file' )
      fullySampledSenseMaps = mri_makeSensitivityMaps( sub_kData, 'mask', mask, 'verbose', verbose );
      save( fullySampledSenseMapsMatFile, 'fullySampledSenseMaps' );
    else
      load( fullySampledSenseMapsMatFile, 'fullySampledSenseMaps' );
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


    fullySampledSakeL1EspiritRecon = sakeRecon2D( sub_kData, 'type', 'espiritL1' );
    if numel( trueRecon ) == 0
      trueRecon = fullySampledSakeL1EspiritRecon;
    else
      printMetrics( logFile, datacase, 1.0, numel(mask), trueRecon, ...
        'trueRecon', datacaseDir, trueRecon );
    end
    printMetrics( logFile, datacase, 1.0, numel(mask), fullySampledSakeL1EspiritRecon, ...
      'sakeL1Espirit', datacaseDir, trueRecon );

    ssqRecon_fullySampled = mri_reconSSQ( sub_kData, 'multiSlice', true);
    printMetrics( logFile, datacase, 1.0, numel(mask), ssqRecon_fullySampled, ...
      'ssqFullySampled', datacaseDir, trueRecon, 'senseMaps', fullySampledSenseMaps );

    parfor sampleFractionIndx = 1 : numel( sampleFractions )
      sampleFraction = sampleFractions( sampleFractionIndx );
      lambda_x = lambda_xs( sampleFractionIndx );
      lambda_s = lambda_ss( sampleFractionIndx );
      lambda_h = lambda_hs( sampleFractionIndx );

      outDir = [ datacaseDir, '/', num2str(sampleFraction,'%3.2f') ];
      if ~exist( outDir, 'dir' ), mkdir( outDir ); end

      [ Ny, Nx, nSlices, nCoils ] = size( sub_kData );
      sMask = [ Ny Nx ];
      if sampleFraction == 1
        vdMask = ones( sMask );
      else
        sampleGoal = round( sampleFraction*prod(sMask) );
        vdMask = vdSampleMask( sMask, vdSigma * [Ny, Nx], sampleGoal );
      end
      kMask = sum( sum( abs( sub_kData ), 3), 4 ) ~= 0;
      vdMask = vdMask & kMask;
      maskFig = figure;  imshowscale( vdMask, 5 );
      saveas( maskFig, [outDir, '/sampleMask.png'] );
      close( maskFig );
      nSamples = sum( vdMask(:) );

      kData = sub_kData;
      for slice=1:nSlices
        for coil=1:nCoils
          kData(:,:,slice,coil) = vdMask .* sub_kData(:,:,slice,coil);
        end
      end

      mccsDir = [ outDir, '/mccs' ];
      [mccsRecon,mccsMaps] = mri_reconMCCS( kData, lambda_x, lambda_s, lambda_h, noiseCoords, ...
        'kcf', kcf, 'res', res, 'maxOuterIter', maxOuterIter, 'noiseCov', noiseCov, ...
        'verbose', verbose, 'doCheckAdjoint', doCheckAdjoint, 'outDir', mccsDir );
      printMetrics( logFile, datacase, sampleFraction, nSamples, mccsRecon, 'mccs', ...
        outDir, trueRecon, 'senseMaps', mccsMaps );
      
      ssqRecon = mri_ssqRecon( kData, 'multiSlice', true );
      printMetrics( logFile, datacase, sampleFraction, nSamples, ssqRecon, ...
        'undersampledSsq', outDir, trueRecon );

      [iSENSEnn_Recon,iSENSEnn_Maps] = mri_iSENSEnn( kData, lambda_iSENSEnn_x, lambda_iSENSEnn_s, ...
        'verbose', verbose, 'doCheckAdjoint', doCheckAdjoint );
      printMetrics( logFile, datacase, sampleFraction, nSamples, iSENSEnn_Recon, 'iSENSEnn', ...
        outDir, trueRecon, 'senseMaps', iSENSEnn_Maps );

      senseLoraksRecon = mri_senseLoraks( kData );
      printMetrics( logFile, datacase, sampleFraction, nSamples, senseLoraksRecon, ...
        'senseLoraks', outDir, trueRecon );

      sakeL1EspiritRecon = sakeRecon2D( kData, 'type', 'espiritL1' );
      printMetrics( logFile, datacase, sampleFraction, nSamples, sakeL1EspiritRecon, ...
        'sakeL1Espirit', outDir, trueRecon );

      [sparseSenseRecon,senseMaps] = mri_csSenseRecon( kData, lambda_csSENSE, 'verbose', verbose );
      printMetrics( logFile, datacase, sampleFraction, nSamples, sparseSenseRecon, 'sparseSense', ...
        outDir, trueRecon, 'senseMaps', senseMaps );

      [nnSenseRecon,nnMaps] = mri_nnSenseRecon( kData, lambda_nnSENSE, 'senseMaps', senseMaps, ...
        'verbose', verbose );
      printMetrics( logFile, datacase, sampleFraction, nSamples, nnSenseRecon, 'nnSense', ...
        outDir, trueRecon, 'senseMaps', nnMaps );

      sakeEspiritRecon = sakeRecon2D( kData, 'type', 'espirit' );
      printMetrics( logFile, datacase, sampleFraction, nSamples, sakeEspiritRecon, ...
        'sakeEspiritRecon', outDir, trueRecon );
    end
  end
end

