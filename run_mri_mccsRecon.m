
function run_mri_mccsRecon
  close all; clear; rng(2);

  datacase = 1;
  doCheckAdjoint = true;
  simCase = 1;
  doSimulation = false;
  lambda_x = 1d-5;
  lambda_s = 1d-2;
  verbose = true;

  sampleFraction = 0.2;
  vdSigma = 75;

  senseMaps = [];
  if doSimulation == true
    [ kData, senseMaps ] = loadSimData( simCase );
  else
    %kData = loadDatacase( datacase );
kData = [];  load kData.mat;
kData = kData( :, :, 4, : );
  end

  kData = kData ./ max( abs( kData(:) ) );

  outDir = [ './out/datacase_', indx2str(datacase,100), '_', num2str(sampleFraction,'%3.2f') ];
  if ~exist( outDir, 'dir' ), mkdir( outDir ); end

  ssqRecon_fullySampled = mri_ssqRecon( kData );
  ssqRecon_fullySampledFig = figure;  imshowscale( abs( ssqRecon_fullySampled ), 5 );
  saveas( ssqRecon_fullySampledFig, [outDir, '/recon_ssqFullySampled.png'] );
  close( ssqRecon_fullySampledFig );

  mask = mri_makeIntensityMask( kData );
  %fullySampledSenseMaps = mri_makeSensitivityMaps( kData, 'mask', mask, 'verbose', verbose );
load fullySampledSenseMaps.mat
  fullySampledSenseRecons = bsxfun( @times, fullySampledSenseMaps, ssqRecon_fullySampled );
  fullySampledSenseReconFig = figure;  showImageCube( abs(fullySampledSenseRecons), 5 );
  titlenice( 'Sense Images - Fully Sampled' );
  saveas( fullySampledSenseReconFig, [outDir, '/recon_fullySampledSenseImages.png'] );
  close( fullySampledSenseReconFig );

  [ Nx, Ny, nSlices, nCoils ] = size( kData );
  sMask = [ Nx Ny ];
  if sampleFraction == 1
    vdMask = ones( sMask );
  else
    vdMask = vdSampleMask( sMask, vdSigma, sampleFraction*prod(sMask) );
  end
  maskFig = figure;  imshowscale( vdMask, 5 )
  saveas( maskFig, [outDir, '/sampleMask.png'] );
  close( maskFig );


  for slice=1:nSlices
    for coil=1:nCoils
      kData(:,:,slice,coil) = vdMask .* kData(:,:,slice,coil);
    end
  end

%   ssqRecon = mri_ssqRecon( kData );
%   undersampledSsqFig = figure;
%   imshowscale( abs( ssqRecon ), 5 );  titlenice( 'Undersampled ssqRecon' );
%   saveas( undersampledSsqFig, [outDir, '/recon_ssq.png'] );
%   close( undersampledSsqFig );
%
% senseMapFile = ['senseMaps_', num2str(sampleFraction), '_.mat'];
% if exist( senseMapFile, 'file' ), load( senseMapFile ); end
%   [csSenseRecon,senseMaps] = mri_csSenseRecon( kData, 100 * lambda_x, 'senseMaps', senseMaps, ...
%     'verbose', verbose );
% if ~exist( senseMapFile, 'file' ), save( senseMapFile, 'senseMaps' ); end
%   csSenseFig = figure; imshowscale( abs(csSenseRecon), 5 );  title('csSense Recon');
%   saveas( csSenseFig, [outDir, '/recon_csSense.png'] );
%   close( csSenseFig );
% 
%   nnSenseRecon = mri_nnSenseRecon( kData, lambda_x*1000, 'senseMaps', senseMaps, ...
%     'verbose', verbose );
%   nnSenseFig = figure; imshowscale( abs(nnSenseRecon), 5 );  title('nnSense Recon');
%   saveas( nnSenseFig, [outDir, '/recon_nnSense.png'] );
%   close( nnSenseFig );

  [iSENSEnn_Recon,iSENSEnn_Maps] = mri_iSENSEnn( kData, lambda_x, lambda_s, ...
    'verbose', verbose, 'doCheckAdjoint', doCheckAdjoint, 'outDir', [ outDir, '/iSENSEnn/' ] );   %#ok<ASGLU>
  iSENSEnn_ReconFig = figure; imshowscale( abs(iSENSEnn_Recon), 5 );  title('mccs Recon');
  saveas( iSENSEnn_ReconFig, [outDir, '/recon_iSENSEnn.png'] );
  close( iSENSEnn_ReconFig );

  [mccsRecon,mccsMaps] = mri_mccsRecon( kData, lambda_x, lambda_s, ...
    'verbose', verbose, 'doCheckAdjoint', doCheckAdjoint, 'outDir', [ outDir, '/mccs/' ] );
  mccsReconFig = figure; imshowscale( abs(mccsRecon), 5 );  title('mccs Recon');
  saveas( mccsReconFig, [outDir, '/recon_mccs.png'] );
  close( mccsReconFig );

  recons = cat( 3, ...
    scaleImg( ssqRecon_fullySampled / max( ssqRecon_fullySampled(:) ), [0 1] ), ...
    scaleImg( ssqRecon / max( ssqRecon(:) ), [0 1] ), ...
    scaleImg( senseRecon / max( senseRecon(:) ), [0 1] ), ...
    scaleImg( mccsRecon / max( mccsRecon(:) ), [0 1] ) );


  figure; showImageCube( senseMaps, 'border', 1, 'borderValue', 'max' );  titlenice( 'sense maps' )
  figure; showImageCube( mccsMaps, 'border', 1, 'borderValue', 'max' );  titlenice( 'mccs maps' );
  figure; showImageCube( recons, 5, 'border', 1, 'borderValue', 'max', 'nImgsPerRow', 4 );
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

