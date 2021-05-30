function IntMatrix = CheckerCalibrate(string)
calibimgset =imageSet(fullfile(string));
imageFileNames = calibimgset.ImageLocation;
[imagePoints, boardSize] = detectCheckerboardPoints(imageFileNames);
squareSizeInMM = 16;
worldPoints = generateCheckerboardPoints(boardSize,squareSizeInMM);
I = readimage(calibimgset,1);
imageSize = [size(I, 1),size(I, 2)];
params = estimateCameraParameters(imagePoints,worldPoints, ...
    'ImageSize',imageSize);
%showReprojectionErrors(params);
%figure;
%showExtrinsics(params);
IntMatrix=params.IntrinsicMatrix;
end