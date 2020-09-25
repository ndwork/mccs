

function out = sakeRecon2D( kData, varargin )
  % out = sakeRecon2D( kData [, 'type', type ] )
  %
  % Inputs:
  % kData is [ nx X ny X nSlices X nc ]
  %
  % Optional Inputs:
  % type - default is espiritL1; alternative is espirit
  %
  % Output - reconstructed volume
  %
  % Written by Nicholas Dwork - Copyright 2019
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  p = inputParser;
  p.addParameter( 'type', 'espirit', @(x) true );
  p.parse( varargin{:} );
  type = p.Results.type;

  origPath = path();
  addpath( genpath( '.' ) );

  nSlices = size( kData, 3 );

  out = cell( 1, 1, nSlices );
  parfor sliceIndx = 1 : nSlices
    out{1,1,sliceIndx} = sakeReconSlice( squeeze( kData(:,:,sliceIndx,:) ), type );   %#ok<PFBNS>
  end
  out = cell2mat( out );

  path( origPath );  % restore original path
end


function out = sakeReconSlice( kData, type )
  % Input: kData is [nx X ny X nc]

  ncalib = 48;
  ksize = [6,6]; % ESPIRiT kernel-window-size
  sakeIter = 100;
  wnthresh = 1.8; % Window-normalized number of singular values to threshold
  eigThresh_im = 0.9; % threshold of eigenvectors in image space
  splitWeight = 0.4;  % reasonable value
  lambda = 0.0025;    % L1-Wavelet threshold
  nIterSplit = 15;    % number of splitting iterations for CS part

  [sx,sy,Nc] = size(kData);

  mask = kData ~= 0;
  kData = kData / max(max(max(abs(ifft2c(kData))))) + eps;
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

  if strcmp( type, 'espirit' )
    [~, out] = cgESPIRiT( kDataC, ESP, nIterCG, 0.01, kDataC*0 );

  elseif strcmp( type, 'espiritL1' )
    XOP = Wavelet('Daubechies_TI',4,6);
    FT = p2DFT( mask, [sx,sy,Nc] );
    tmp = zeros( size( kDataC, 1 ), size( kDataC, 2 ), 2 );
    out = cgL1ESPIRiT( kDataC, tmp, FT, ESP, nIterCG, XOP, lambda, splitWeight, nIterSplit );

  end

  out = out(:,:,2);
end
