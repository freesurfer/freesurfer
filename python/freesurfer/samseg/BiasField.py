import numpy as np
import math

eps = np.finfo( float ).eps


class BiasField:
    def __init__(self, imageSize, smoothingKernelSize):
        self.fullBasisFunctions = self.getBiasFieldBasisFunctions(imageSize, smoothingKernelSize)
        self.basisFunctions = self.fullBasisFunctions.copy()
        self.coefficients = None

    def backprojectKroneckerProductBasisFunctions(self, kroneckerProductBasisFunctions, coefficients):
        numberOfDimensions = len(kroneckerProductBasisFunctions)
        Ms = np.zeros(numberOfDimensions, dtype=np.uint32)  # Number of basis functions in each dimension
        Ns = np.zeros(numberOfDimensions, dtype=np.uint32)  # Number of basis functions in each dimension
        transposedKroneckerProductBasisFunctions = []
        for dimensionNumber in range(numberOfDimensions):
            Ms[dimensionNumber] = kroneckerProductBasisFunctions[dimensionNumber].shape[1]
            Ns[dimensionNumber] = kroneckerProductBasisFunctions[dimensionNumber].shape[0]
            transposedKroneckerProductBasisFunctions.append(kroneckerProductBasisFunctions[dimensionNumber].T)
        y = self.projectKroneckerProductBasisFunctions(transposedKroneckerProductBasisFunctions,
                                                       coefficients.reshape(Ms, order='F') )
        Y = y.reshape(Ns, order='F')
        return Y

    def projectKroneckerProductBasisFunctions(self, kroneckerProductBasisFunctions, T):
        #
        # Compute
        #   c = W' * t
        # where
        #   W = W{ numberOfDimensions } \kron W{ numberOfDimensions-1 } \kron ... W{ 1 }
        # and
        #   t = T( : )
        numberOfDimensions = len(kroneckerProductBasisFunctions)
        currentSizeOfT = list(T.shape)
        for dimensionNumber in range(numberOfDimensions):
            # Reshape into 2-D, do the work in the first dimension, and shape into N-D
            T = T.reshape((currentSizeOfT[0], -1), order='F')
            T = ( kroneckerProductBasisFunctions[dimensionNumber] ).T @ T
            currentSizeOfT[0] = kroneckerProductBasisFunctions[dimensionNumber].shape[1]
            T = T.reshape(currentSizeOfT, order='F')
            # Shift dimension
            currentSizeOfT = currentSizeOfT[1:] + [currentSizeOfT[0]]
            T = np.rollaxis(T, 0, 3)
        # Return result as vector
        coefficients = T.flatten(order='F')
        return coefficients

    def computePrecisionOfKroneckerProductBasisFunctions(self, kroneckerProductBasisFunctions, B):
        #
        # Compute
        #   H = W' * diag( B ) * W
        # where
        #   W = W{ numberOfDimensions } \kron W{ numberOfDimensions-1 } \kron ... W{ 1 }
        # and B is a weight matrix
        numberOfDimensions = len( kroneckerProductBasisFunctions )

        # Compute a new set of basis functions (point-wise product of each combination of pairs) so that we can
        # easily compute a mangled version of the result
        Ms = np.zeros( numberOfDimensions, dtype=np.uint32) # Number of basis functions in each dimension
        hessianKroneckerProductBasisFunctions = {}
        for dimensionNumber in range(numberOfDimensions):
            M = kroneckerProductBasisFunctions[dimensionNumber].shape[1]
            A = kroneckerProductBasisFunctions[dimensionNumber]
            hessianKroneckerProductBasisFunctions[dimensionNumber] = np.kron( np.ones( (1, M )), A ) * np.kron( A, np.ones( (1, M) ) )
            Ms[dimensionNumber] = M
        result = self.projectKroneckerProductBasisFunctions( hessianKroneckerProductBasisFunctions, B )
        new_shape = list(np.kron( Ms, [ 1, 1 ] ))
        new_shape.reverse()
        result = result.reshape(new_shape)
        permutationIndices = np.hstack((2 * np.r_[: numberOfDimensions ], 2 * np.r_[: numberOfDimensions ] +1))
        result = np.transpose(result, permutationIndices)
        precisionMatrix = result.reshape( ( np.prod( Ms ), np.prod( Ms ) ) )
        return precisionMatrix

    def getBiasFieldBasisFunctions(self, imageSize, smoothingKernelSize):
        # Our bias model is a linear combination of a set of basis functions. We are using so-called
        # "DCT-II" basis functions, i.e., the lowest few frequency components of the Discrete Cosine
        # Transform.
        biasFieldBasisFunctions = []
        for dimensionNumber in range(3):
            N = imageSize[dimensionNumber]
            delta = smoothingKernelSize[dimensionNumber]
            M = math.ceil(N / delta) + 1
            Nvirtual = (M - 1) * delta
            js = [(index + 0.5) * math.pi / Nvirtual for index in range(N)]
            scaling = [math.sqrt(2 / Nvirtual)] * M
            scaling[0] /= math.sqrt(2)
            A = np.array([[math.cos(freq * m) * scaling[m] for m in range(M)] for freq in js])
            biasFieldBasisFunctions.append(A)

        return biasFieldBasisFunctions

    def getBiasFields(self, mask=None):
        #
        numberOfContrasts = self.coefficients.shape[-1]
        imageSize = tuple([functions.shape[0] for functions in self.basisFunctions])
        biasFields = np.zeros(imageSize + (numberOfContrasts,), order='F')
        for contrastNumber in range(numberOfContrasts):
            biasField = self.backprojectKroneckerProductBasisFunctions(
                self.basisFunctions, self.coefficients[:, contrastNumber])
            if mask is not None:
                biasField *= mask
            biasFields[:, :, :, contrastNumber] = biasField

        return biasFields

    def fitBiasFieldParameters(self, imageBuffers, gaussianPosteriors, means, variances, mask):

        # Bias field correction: implements Eq. 8 in the paper
        #    Van Leemput, "Automated Model-based Bias Field Correction of MR Images of the Brain", IEEE TMI 1999

        #
        numberOfGaussians = means.shape[0]
        numberOfContrasts = means.shape[1]
        numberOfBasisFunctions = [functions.shape[1] for functions in self.basisFunctions]
        numberOf3DBasisFunctions = np.prod(numberOfBasisFunctions)

        # Set up the linear system lhs * x = rhs
        precisions = np.zeros_like(variances)
        for gaussianNumber in range(numberOfGaussians):
            precisions[gaussianNumber, :, :] = np.linalg.inv(variances[gaussianNumber, :, :]).reshape(
                (1, numberOfContrasts, numberOfContrasts))

        lhs = np.zeros((numberOf3DBasisFunctions * numberOfContrasts,
                        numberOf3DBasisFunctions * numberOfContrasts))  # left-hand side of linear system
        rhs = np.zeros((numberOf3DBasisFunctions * numberOfContrasts, 1))  # right-hand side of linear system
        weightsImageBuffer = np.zeros(mask.shape)
        tmpImageBuffer = np.zeros(mask.shape)
        for contrastNumber1 in range(numberOfContrasts):
            # logger.debug('third time contrastNumber=%d', contrastNumber)
            contrast1Indices = np.arange(0, numberOf3DBasisFunctions) + \
                               contrastNumber1 * numberOf3DBasisFunctions

            tmp = np.zeros(gaussianPosteriors.shape[0])
            for contrastNumber2 in range(numberOfContrasts):
                contrast2Indices = np.arange(0, numberOf3DBasisFunctions) + \
                                   contrastNumber2 * numberOf3DBasisFunctions

                classSpecificWeights = gaussianPosteriors * precisions[:, contrastNumber1, contrastNumber2]
                weights = np.sum(classSpecificWeights, 1)

                # Build up stuff needed for rhs
                predicted = np.sum(classSpecificWeights * means[:, contrastNumber2], 1) / (weights + eps)
                residue = imageBuffers[mask, contrastNumber2] - predicted
                tmp += weights * residue

                # Fill in submatrix of lhs
                weightsImageBuffer[mask] = weights
                lhs[np.ix_(contrast1Indices, contrast2Indices)] \
                    = self.computePrecisionOfKroneckerProductBasisFunctions(self.basisFunctions,
                                                                       weightsImageBuffer)

            tmpImageBuffer[mask] = tmp
            rhs[contrast1Indices] = self.projectKroneckerProductBasisFunctions(self.basisFunctions,
                                                                          tmpImageBuffer).reshape(-1, 1)

        # Solve the linear system x = lhs \ rhs
        solution = np.linalg.solve(lhs, rhs)

        #
        self.coefficients = solution.reshape((numberOfContrasts, numberOf3DBasisFunctions)).transpose()

    def setBiasFieldCoefficients(self, coefficients):
        self.coefficients = coefficients

    def downSampleBasisFunctions(self, downSamplingFactors):
        self.basisFunctions = [np.array(biasFieldBasisFunction[::downSamplingFactor])
                                        for biasFieldBasisFunction, downSamplingFactor in
                                        zip(self.fullBasisFunctions, downSamplingFactors)]
