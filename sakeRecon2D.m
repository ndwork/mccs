

function out = sakeRecon2D( kData )
  % Inputs:
  % kData is [ nx X ny X nSlices X nc ]

  p = path();
  addpath( genpath( '.' ) );
  
  nSlices = size( kData, 3 );
  
  out = cell( 1, 1, nSlices );
  %parfor sliceIndx = 1 : nSlices
for sliceIndx = 1 : nSlices
    out{1,1,sliceIndx} = sakeReconSlice( squeeze( kData(:,:,sliceIndx,:) ) );   %#ok<PFBNS>
  end
  out = cell2mat( out );

  path( p );  % restore original path
end


function out = sakeReconSlice( kData )
  % Inputs:
  % kData is [nx X ny X nc]

  ncalib = 48;
  ksize = [6,6]; % ESPIRiT kernel-window-size
  sakeIter = 100;
  wnthresh = 1.8; % Window-normalized number of singular values to threshold
  eigThresh_im = 0.9; % threshold of eigenvectors in image space

  [sx,sy,Nc] = size(kData);

  mask = kData ~= 0;
  kData = kData/max(max(max(abs(ifft2c(kData))))) + eps;
  kDataC = kData .* mask;

  calibc = crop( kDataC, [ncalib,ncalib,Nc] );
  calib = SAKE( calibc, [ksize], wnthresh, sakeIter, 0 );


  %% Singular values of the calibration matrix and ESPIRiT Maps after SAKE
  % Sake now shows a null space and improved Maps. 

  [k,S] = dat2Kernel(calib,ksize);   %#ok<ASGLU>
  [M,W] = kernelEig(k(:,:,:,1:floor(wnthresh*prod(ksize))),[sx,sy]);


  %% Compute Soft-SENSE ESPIRiT Maps 
  % crop sensitivity maps according to eigenvalues==1. Note that we have to
  % use 2 sets of maps. Here we weight the 2 maps with the eigen-values

  maps = M(:,:,:,end-1:end);

  % Weight the eigenvectors with soft-senses eigen-values
  weights = W(:,:,end-1:end) ;
  weights = (weights - eigThresh_im)./(1-eigThresh_im).* (W(:,:,end-1:end) > eigThresh_im);
  weights = -cos(pi*weights)/2 + 1/2;

  % create and ESPIRiT operator
  ESP = ESPIRiT( maps, weights );
  nIterCG = 15; 

  [~, out] = cgESPIRiT(kDataC,ESP, nIterCG, 0.01,kDataC*0);
  out = out(:,:,2);

end
