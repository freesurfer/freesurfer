lhs = zeros(nrOfChannels*nrOfUsedBasisFuncs);
rhs = zeros(nrOfChannels*nrOfUsedBasisFuncs, 1);
orthoSampleBasisFunctions = zeros(DIM,nrOfUsedBasisFuncs);
orthoSampleBasisFunctions(maskIndices,:) = ...
reshape(ortho);


precisions = zeros(size(sigmas));
  
for class = 1 : sum(numberOfGaussiansPerClass)
    precisions(:,:,class) = inv(sigmas(:,:,class));
end

predictedImage = zeros( DIM, 2*numberOfInputImages );
weightImage = zeros( DIM, 2*numberOfInputImages  );
dataImage = zeros( DIM, 2*numberOfInputImages );
for channel1 = 1:numberOfInputImages
    for channel2 = 1:numberOfInputImages
        tmp = posteriors ./ repmat( precisions(channel1,channel2,:)', [ size( posteriors, 1 ) 1 ] );
        predictedImage(maskIndices,(channel1-1) + channel2) = sum( tmp .* repmat( means(channel1,:), [ size( posteriors, 1 ) 1 ] ), 2 ) ./ ( sum(tmp,2) + eps );
        weightImage(maskIndices,(channel1-1) + channel2) = sum(tmp, 2 );
        dataImage(maskIndices, (channel1-1) + channel2) = data(:,channel2);
    end
end
predictedImage = reshape(predictedImage,[prod(DIM) 2*numberOfInputImages]);
weightImage = reshape(weightImage,[prod(DIM) 2*numberOfInputImages]);
dataImage = reshape(dataImage,[prod(DIM) 2*numberOfInputImages]);

residueImage = dataImage - predictedImage; % The bias field is simply a (weighted) smoothed version of this

for channel = 1:numberOfInputImages
    rhs((channel-1)*numberOfBasisFunctions + [1:numberOfBasisFunctions]) = orthoSampleBasisFunctions'*sum(weightImage(:,channel:2*channel).*...
                                                                            residueImage(:,channel:2*channel),2);
    for channel2 = 1:numberOfInputImages
        lhs((channel-1)*numberOfBasisFunctions + [1:numberOfBasisFunctions],...
            (channel2-1)*numberOfBasisFunctions + [1:numberOfBasisFunctions]) = ...
            orthoSampledBasisFunction'*(repmat(weightImage(:,(channel-1) + channel1),[1 numberOfBasisFunctions]).*orthoSampleBasisFuncs);
    end
end

solution = lhs\rhs;




for row=1:nrOfChannels
    tmp = squeeze(precisions(row,:,1:end));
    if nrOfChannels==1
        tmp = tmp';
    end
    rhs((row-1)*nrOfUsedBasisFuncs+[1:nrOfUsedBasisFuncs]) = ...
            orthoSampleBasisFuncs(:,1:nrOfUsedBasisFuncs)' * ...
            sum((squeeze(weights(:,row,:)) .* sampleData - ...
            sampleClassification(:,1:end-1) * ...
            (tmp .* means(:,1:end-1))'), 2);
    for col=1:nrOfChannels
        if (col>=row)
            lhs((row-1)*nrOfUsedBasisFuncs+[1:nrOfUsedBasisFuncs], ...
                (col-1)*nrOfUsedBasisFuncs+[1:nrOfUsedBasisFuncs]) = ...
                orthoSampleBasisFuncs(:,1:nrOfUsedBasisFuncs)' * ...
                ((weights(:,row,col) * ones(1, nrOfUsedBasisFuncs)) .* ...
                sampleBasisFuncs(:,1:nrOfUsedBasisFuncs));
        else
            lhs((row-1)*nrOfUsedBasisFuncs+[1:nrOfUsedBasisFuncs], ...
                (col-1)*nrOfUsedBasisFuncs+[1:nrOfUsedBasisFuncs]) = ...
                lhs((col-1)*nrOfUsedBasisFuncs+[1:nrOfUsedBasisFuncs], ...
                (row-1)*nrOfUsedBasisFuncs+[1:nrOfUsedBasisFuncs]);
        end
    end
end
solution = lhs\rhs;
solution = reshape(solution, [nrOfUsedBasisFuncs nrOfChannels]);
for channel=1:nrOfChannels
    biasCoeff(1:nrOfUsedBasisFuncs,channel) = solution(:,channel);
    biasCoeff(nrOfUsedBasisFuncs+1:nrOfBasisFuncs,channel) = 0;
end




            
    precisions = zeros(size(sigmas));
  
    for class = 1 : sum(numberOfGaussiansPerClass)
        precisions(:,:,class) = inv(sigmas(:,:,class));
    end

    predictedImage = zeros( DIM, 2*numberOfInputImages );
    weightImage = zeros( DIM, 2*numberOfInputImages  );
    dataImage = zeros( DIM, 2*numberOfInputImages );
    for channel1 = 1:numberOfInputImages
        for channel2 = 1:numberOfInputImages
        
            tmp = posteriors ./ repmat( precisions(channel1,channel2,:)', [ size( posteriors, 1 ) 1 ] );
            predictedImage(maskIndices,(channel1-1) + channel2) = sum( tmp .* repmat( means(channel1,:), [ size( posteriors, 1 ) 1 ] ), 2 ) ./ ( sum(tmp,2) + eps );
            weightImage(maskIndices,(channel1-1) + channel2) = sum(tmp, 2 );
            dataImage(maskIndices, (channel1-1) + channel2) = data(:,channel2);
        end
    end
    
    residueImage = dataImage - predictedImage; % The bias field is simply a (weighted) smoothed version of this
    if ( EMIterationNumber == 1 )
        [ estimatedBiasField, coefficients ] = smoothWithSeparableBasisFunctionsWithWeightsNImages( residueImage, weightImage, numberOfBasisFunctions, maximumNumberOfBasisFunctions, nima, numberOfInputImages );
    else
        [ estimatedBiasField, coefficients ] = smoothWithSeparableBasisFunctionsWithWeightsNImages( residueImage, weightImage, numberOfBasisFunctions, nima, numberOfInputImages ); % This is just a faster version of the above, but needs the above to have been
                                    % called at least once...
    end
            
    estimatedBiasFields(:,:,:,nima) = estimatedBiasField;
    biasCorrectedData(:,nima) = data(:,nima) - estimatedBiasField( maskIndices); % Correct the data for the estimated bias field
        
    if(numberOfBasisFunctions == maximumNumberOfBasisFunctions)
          saveCoefficients(:,:,:,nima) = coefficients;
    end

