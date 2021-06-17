
function run_reconstructions
  close all; clear; rng(1);

  datacases = [ 1, 0, 11, 2, 6, 5, 9, 3, 7, 8 ];
  % Note, datacase 0 means simulation

  outDir = './out/';
  doCheckAdjoint = false;
  simCase = 2;
  vdSigma = 0.3;
  maxOuterIter = 50;

  sampleFractions = [ 0.10 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.75 ];
  verbose = true;

  if ~exist( 'sub_kData', 'dir' ), mkdir( 'sub_kData' ); end
  if ~exist( 'fullySampledSenseMaps', 'dir' ), mkdir( 'fullySampledSenseMaps' ); end

  for datacase = datacases
    disp([ 'Working on datacase ', num2str( datacase ) ]);

    delete( gcp );
    parpool( 16 );

    reconstructDatacase( datacase, sampleFractions, outDir, ...
      doCheckAdjoint, simCase, vdSigma, maxOuterIter, verbose );
  end
end


function reconstructDatacase( datacase, sampleFractions, outDir, ...
  doCheckAdjoint, simCase, vdSigma, maxOuterIter, verbose )

  dataMatFile = [ 'sub_kData/sub_kData_', indx2str(datacase,99), '.mat' ];
  if exist( dataMatFile, 'file' )
    loadData = false;
    [ ~, noiseCoords, kcfs, res, lambda_xs, lambda_ss, lambda_hs, lambda_iSENSEnn_xs, ...
      lambda_iSENSEnn_ss, lambda_sparseSENSEs, lambda_nnSENSEs, lambda_espiritL1s, ...
      trueRecon ] = loadDatacase( datacase, loadData, 'simCase', simCase );
    load( dataMatFile, 'sub_kData' );
  else
    loadData = true;
    [ sub_kData, noiseCoords, kcfs, res, lambda_xs, lambda_ss, lambda_hs, lambda_iSENSEnn_xs, ...
      lambda_iSENSEnn_ss, lambda_sparseSENSEs, lambda_nnSENSEs, lambda_espiritL1s, ...
      trueRecon ] = loadDatacase( datacase, loadData, 'simCase', simCase );
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

  datacaseDir = [ outDir, '/case_', indx2str(datacase,100) ];
  if ~exist( datacaseDir, 'dir' ), mkdir( datacaseDir ); end
  logFile = [ datacaseDir, '/metrics.log' ];
  logID = fopen( logFile, 'w' );
  fprintf( logID, [ 'datacase, sampleFraction, nSamples, type, mdm, mdm-2, niqe, ', ...
    'piqe, mse, mae, ssim, correlation, angleErr, mi \n' ] );
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

  [ fullySampledSakeL1EspiritRecon, fullySampledSakeL1EspiritMaps ] = loadSavedResult( ...
    datacaseDir, 'sakeL1Espirit' );
  if numel( fullySampledSakeL1EspiritMaps ) == 0
    [ fullySampledSakeL1EspiritRecon, fullySampledSakeL1EspiritMaps ] = sakeRecon2D( ...
      sub_kData, 'type', 'espiritL1' );
  end
  if numel( trueRecon ) == 0
    trueRecon = fullySampledSakeL1EspiritRecon;
  end
  printMetrics( logFile, datacase, 1.0, numel(mask), fullySampledSakeL1EspiritRecon, ...
    'sakeL1Espirit', datacaseDir, trueRecon, 'senseMaps', fullySampledSakeL1EspiritMaps );
  clear fullySampledSakeL1EspiritRecon fullySampledSakeL1EspiritMaps

  fullySampledCoilRecons = mri_fftRecon( sub_kData, 'multiSlice', true );
  roemer_fullySampled = mri_reconRoemer( fullySampledCoilRecons );
  printMetrics( logFile, datacase, 1.0, numel(mask), roemer_fullySampled, ...
    'roemerFullySampled', datacaseDir, trueRecon, 'senseMaps', fullySampledSenseMaps );
  clear fullySampledCoilRecons roemer_fullySampled

  ssqRecon_fullySampled = mri_reconSSQ( sub_kData, 'multiSlice', true );
  printMetrics( logFile, datacase, 1.0, numel(mask), ssqRecon_fullySampled, ...
    'ssqFullySampled', datacaseDir, trueRecon, 'senseMaps', fullySampledSenseMaps );
  clear ssqRecon_fullySampled

  parfor sampleFractionIndx = 1 : numel( sampleFractions )
    sampleFraction = sampleFractions( sampleFractionIndx );
    disp(['Working on sample fraction ', num2str( sampleFraction ) ] );

    reconstructSampleFraction( logFile, datacase, datacaseDir, sampleFraction, sub_kData, kcfs, res, ...
      lambda_xs, lambda_ss, lambda_hs, lambda_iSENSEnn_xs, lambda_iSENSEnn_ss, lambda_sparseSENSEs, ...
      lambda_nnSENSEs, lambda_espiritL1s, trueRecon, noiseCov, noiseCoords, maxOuterIter, vdSigma, doCheckAdjoint, ...
      verbose );
  end

end


function reconstructSampleFraction( logFile, datacase, datacaseDir, sampleFraction, sub_kData, kcfs, res, lambda_xs, lambda_ss, lambda_hs, ...
  lambda_iSENSEnn_xs, lambda_iSENSEnn_ss, lambda_sparseSENSEs, lambda_nnSENSEs, ...
  lambda_espiritL1s, trueRecon, noiseCov, noiseCoords, maxOuterIter, vdSigma, doCheckAdjoint, verbose )

  outSampleDir = [ datacaseDir, '/', num2str(sampleFraction,'%3.2f') ];
  if ~exist( outSampleDir, 'dir' ), mkdir( outSampleDir ); end

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
  saveas( maskFig, [ outSampleDir, '/sampleMask.png' ] );
  close( maskFig );
  nSamples = sum( vdMask(:) );

  kData = sub_kData;
  for slice=1:nSlices
    for coil=1:nCoils
      kData(:,:,slice,coil) = vdMask .* sub_kData(:,:,slice,coil);
    end
  end

  roemerCoilRecons = mri_fftRecon( kData, 'multiSlice', true );
  roemerRecon = mri_reconRoemer( roemerCoilRecons );
  printMetrics( logFile, datacase, sampleFraction, nSamples, roemerRecon, ...
    'roemer', outSampleDir, trueRecon );
  clear roemerCoilRecons roemerRecon

  ssqRecon = mri_reconSSQ( kData, 'multiSlice', true );
  printMetrics( logFile, datacase, sampleFraction, nSamples, ssqRecon, ...
    'undersampledSsq', outSampleDir, trueRecon );
  clear ssqRecon

  disp( 'Running iSENSEnn' );
  for indxLambdaSs_iSENSE = 1 : numel( lambda_iSENSEnn_ss )
    lambda_iSENSEnn_s = lambda_iSENSEnn_ss( indxLambdaSs_iSENSE );
    for indxLambdaXs_iSENSE = 1 : numel( lambda_iSENSEnn_xs )
      lambda_iSENSEnn_x = lambda_iSENSEnn_xs( indxLambdaXs_iSENSE );

      iSENSE_Dir = [ outSampleDir, '/iSENSE/Ls_', num2str( lambda_iSENSEnn_s, '%3.2e' ), ...
        '/Lx_', num2str( lambda_iSENSEnn_x, '%3.2e' ) ];
      if ~exist( iSENSE_Dir, 'dir' ), mkdir( iSENSE_Dir ); end

      [iSENSEnn_Recon,iSENSEnn_Maps] = loadSavedResult( iSENSE_Dir, 'iSENSEnn' );
      if numel( iSENSEnn_Maps ) == 0
        [iSENSEnn_Recon,iSENSEnn_Maps] = mri_iSENSEnn( kData, lambda_iSENSEnn_x, lambda_iSENSEnn_s, ...
          'verbose', verbose, 'doCheckAdjoint', doCheckAdjoint );
      end
      printMetrics( logFile, datacase, sampleFraction, nSamples, iSENSEnn_Recon, 'iSENSEnn', ...
        iSENSE_Dir, trueRecon, 'senseMaps', iSENSEnn_Maps, 'optLogNames', { 'lambda_s', 'lambda_x' }, ...
        'optLogValues', [ lambda_iSENSEnn_s, lambda_iSENSEnn_x ] );
    end
  end
  clear iSENSEnn_Recon iSENSEnn_Maps

  disp( 'Running senseLoraks' );
  senseLoraksRecon = loadSavedResult( outSampleDir, 'senseLoraks' );
  if numel( senseLoraksRecon ) == 0
    senseLoraksRecon = mri_senseLoraks( kData );
  end
  printMetrics( logFile, datacase, sampleFraction, nSamples, senseLoraksRecon, ...
    'senseLoraks', outSampleDir, trueRecon );
  clear senseLoraksRecon

  disp( 'Running sakeEspirit' );
  [ sakeEspiritRecon, sakeEspiritMaps ]  = loadSavedResult( outSampleDir, 'sakeEspiritRecon' );
  if numel( sakeEspiritMaps ) == 0
    [ sakeEspiritRecon, sakeEspiritMaps ] = sakeRecon2D( kData, 'type', 'espirit' );
  end
  printMetrics( logFile, datacase, sampleFraction, nSamples, sakeEspiritRecon, ...
    'sakeEspiritRecon', outSampleDir, trueRecon, 'senseMaps', sakeEspiritMaps );
  clear sakeEspiritRecon sakeEspiritMaps

  disp( 'Running sakeL1Espirit' );
  for indxLambda_espiritL1 = 1 : numel( lambda_espiritL1s )
    lambda_espiritL1 = lambda_espiritL1s( indxLambda_espiritL1 );

    espiritL1Dir = [ outSampleDir, '/espiritL1/L_', num2str( lambda_espiritL1, '%3.2e' ) ];
    if ~exist( espiritL1Dir, 'dir' ), mkdir( espiritL1Dir ); end

    [ sakeL1EspiritRecon, sakeL1EspiritMaps ] = loadSavedResult( espiritL1Dir, 'sakeL1Espirit' );
    if numel( sakeL1EspiritMaps ) == 0
      [ sakeL1EspiritRecon, sakeL1EspiritMaps ] = sakeRecon2D( ...
        kData, 'type', 'espiritL1', 'lambda', lambda_espiritL1 );
    end
    printMetrics( logFile, datacase, sampleFraction, nSamples, sakeL1EspiritRecon, ...
      'sakeL1Espirit', espiritL1Dir, trueRecon, 'senseMaps', sakeL1EspiritMaps, ...
      'optLogNames', {'lambda_espiritL1'}, 'optLogValues', lambda_espiritL1 );
  end
  clear sakeL1EspiritRecon sakeL1EspiritMaps

  disp( 'Running sparseSense' );
  for indxLambda_sparseSENSE = 1 : numel( lambda_sparseSENSEs )
    lambda_sparseSENSE = lambda_sparseSENSEs( indxLambda_sparseSENSE );

    sparseSENSE_Dir = [ outSampleDir, '/sparseSENSE/L_', num2str( lambda_sparseSENSE, '%3.2e' ) ];
    if ~exist( sparseSENSE_Dir, 'dir' ), mkdir( sparseSENSE_Dir ); end

    [sparseSenseRecon,sparseSenseMaps] = loadSavedResult( sparseSENSE_Dir, 'sparseSense' );
    if numel( sparseSenseMaps ) == 0
      [sparseSenseRecon,sparseSenseMaps] = mri_sparseSENSERecon( kData, lambda_sparseSENSE, ...
        'verbose', verbose );
    end
    printMetrics( logFile, datacase, sampleFraction, nSamples, sparseSenseRecon, 'sparseSense', ...
      sparseSENSE_Dir, trueRecon, 'senseMaps', sparseSenseMaps );
  end
  clear sparseSenseRecon

  disp( 'Running nnSENSE' );
  for indxLambda_nnSENSE = 1 : numel( lambda_nnSENSEs )
    lambda_nnSENSE = lambda_nnSENSEs( indxLambda_nnSENSE );

    nnSENSE_Dir = [ outSampleDir, '/nnSENSE/L_', num2str( lambda_nnSENSE, '%3.2e' ) ];
    if ~exist( nnSENSE_Dir, 'dir' ), mkdir( nnSENSE_Dir ); end

    [ nnSenseRecon, nnSenseMaps ] = loadSavedResult( nnSENSE_Dir, 'nnSense' );
    if numel( nnSenseMaps ) == 0
      [nnSenseRecon,nnSenseMaps] = mri_nnSenseRecon( kData, lambda_nnSENSE, ...
        'senseMaps', sparseSenseMaps, 'verbose', verbose );
    end
    printMetrics( logFile, datacase, sampleFraction, nSamples, nnSenseRecon, 'nnSense', ...
      nnSENSE_Dir, trueRecon, 'senseMaps', nnSenseMaps, 'optLogNames', {'lambda_nnSENSE'}, ...
      'optLogValues', lambda_nnSENSE );
    clear nnSenseRecon nnSenseMaps 
  end
  clear sparseSenseMaps

  disp( 'Running MCCS' );
  runReconsMCCS( trueRecon, datacase, sampleFraction, nSamples, kData, kcfs, lambda_ss, lambda_hs, lambda_xs, outSampleDir, noiseCoords, ...
    res, maxOuterIter, noiseCov, logFile, verbose, doCheckAdjoint );

end


function runReconsMCCS( trueRecon, datacase, sampleFraction, nSamples, kData, kcfs, lambda_ss, lambda_hs, lambda_xs, outSampleDir, noiseCoords, ...
  res, maxOuterIter, noiseCov, logFile, verbose, doCheckAdjoint )
  for kcfIndx = 1 : numel( kcfs )
    kcf = kcfs( kcfIndx );

    for indxLambdaS = 1 : numel( lambda_ss )
      lambda_s = lambda_ss( indxLambdaS );
      for indxLambdaH = 1 : numel( lambda_hs )
        lambda_h = lambda_hs( indxLambdaH );
        for indxLambdaX = 1 : numel( lambda_xs )
          lambda_x = lambda_xs( indxLambdaX );

          mccsDir = [ outSampleDir, '/mccs/kcf_', num2str(kcf), ...
            '/Ls_', num2str( lambda_s, '%3.2e' ), ...
            '/Lh_', num2str( lambda_h, '%3.2e' ), ...
            '/Lx_', num2str( lambda_x, '%3.2e' ) ];
          if ~exist( mccsDir, 'dir' ), mkdir( mccsDir ); end

          [ mccsRecon, mccsMaps ] = loadSavedResult( mccsDir, 'mccs' );
          if numel( mccsMaps ) == 0
            [mccsRecon,mccsMaps] = mri_reconMCCS( kData, lambda_x, lambda_s, lambda_h, noiseCoords, ...
              'kcf', kcf, 'res', res, 'maxOuterIter', maxOuterIter, 'noiseCov', noiseCov, ...
              'verbose', verbose, 'doCheckAdjoint', doCheckAdjoint, 'outDir', mccsDir );
          end
          optLogNames = { 'lambda_s', 'lambda_h', 'lambda_x' };
          optLogValues = [ lambda_s, lambda_h, lambda_x ];
          printMetrics( logFile, datacase, sampleFraction, nSamples, mccsRecon, 'mccs', ...
            mccsDir, trueRecon, 'senseMaps', mccsMaps, 'optLogNames', optLogNames, ...
            'optLogValues', optLogValues );
          clear mccsRecon mccsMaps
        end
      end
    end
  end
end

