
outputFileName = 'Onion02.tif'
for K=1:length(imOut3D(1, 1, :))
   imwrite(uint8(imOut3D(:, :, K)), outputFileName, 'WriteMode', 'append');
end