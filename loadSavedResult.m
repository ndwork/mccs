
function [ recon, senseMaps ] = loadSavedResult( outDir, imgTitle )

  recon = [];
  senseMaps = [];

  saveFile = [ outDir, '/mat_recon_', imgTitle, '.mat' ];
  if exist( saveFile, 'file' ) == 0, return; end

  savedVariables = who( '-file', saveFile );

  if ismember( 'senseMaps', savedVariables )
    load( saveFile, 'recon', 'senseMaps' );
  else
    load( saveFile, 'recon' );
  end

end

