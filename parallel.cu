#include <iostream>
#include <curand_kernel.h>
#include <chrono>

using namespace std;

void sieveOfEratosthenes(int *sieveInput, int max) {
    for(int i = 0; i < max-1; ++i) {
        sieveInput[i] = i+2;
    }
    for(int j= 0; j < max-1; ++j) {
        // only begin sieve round if starting point is not 0
        if(sieveInput[j]) {
            for(int k = j+sieveInput[j]; k < max-1; k += sieveInput[j]) {
                if(sieveInput[k]) {
                    sieveInput[k] = 0;
                }
            }
        }
    }
}

// remove zeroes from sieve, will only contain prime numbers
void compressSieveOfEratosthenes(int *uncompSieve, int *compSieve, int uncompressedSieveSize) {
    int counter = 0;
    for(int i = 0; i < uncompressedSieveSize; ++i) {
        if(uncompSieve[i] != 0) {
            compSieve[counter] = uncompSieve[i];
            counter++;
        }
    }
}

void generateSieve(int *uncompSieve, int maxPrime) {
    sieveOfEratosthenes(uncompSieve, maxPrime);
}

void generateCompressedSieve(int *uncompSieve, int *compSieve, int maxPrime) {
    compressSieveOfEratosthenes(uncompSieve, compSieve, maxPrime);
}

// compressed size will be used in main and kernel functions
int getCompressedSize(int *sieve, int maxPrime) {
    int compressedSize = 0;
    for(int i = 0; i < maxPrime; ++i) {
        if(sieve[i] != 0) {
            compressedSize++;
        }
    }

    return compressedSize - 1;
}

__device__ int getRandomNumber(int idx, int min, int max, curandState state) {
    curand_init(clock64(), idx, 0, &state);
    int difference = max - min + 1;
    int randomNum = (float)curand_uniform(&state)*difference + min;

    return randomNum;
}

// generate numerical representation of word
__device__ int generateWordNum(int size, int idx, curandState state) {
    curand_init(clock64(), idx, 0, &state);
    int wordNum = 0;
    int counter = 0;
    for(int j = 0; j < size; ++j) {
        int randomNum = getRandomNumber(idx, 1, 26, state);
        wordNum = wordNum + randomNum * pow(27, counter);
        counter++;
    }

    return wordNum;
}

__device__ int extendedEuclidGcd(long int first, long int second, long int *x, long int *y) {
    if(!first) {
        *x = 0;
        *y = 1;
        return second;
    }

    long int x1;
    long int y1;
    int gcdResult = extendedEuclidGcd(second%first, first, &x1, &y1);

    *x = y1 - (floorf(second/first)) * x1;
    *y = x1;

    return gcdResult;
}

__device__ long int euclidGcd(long int first, long int second) {
    if(!second) {
        return first;
    } else {
        long int remainder = first%second;
        return euclidGcd(second, remainder);
    }
}

// e is an odd number that is relatively prime to Phi-n
__device__ int generateE(int idx, curandState state, long int phiN) {
    int e = 2;
    while(e == 2 || e%2 == 0 || euclidGcd(e, phiN) != 1) {
        e = getRandomNumber(idx, 100, 500, state);
    }

    return e;
}

// calculates modulus of exponentiation without overflow
__device__ int exponentiationRemainder(long int a, long int b, long int c) {
    int result = 1;

    while(b > 0) {
        if(b%2 == 1) {
            result = (result*a)%c;
        }
        b = floorf(b/2);
        a = (a*a)%c;
    }

    return result;
}

__device__ int encrypt(long int num, long int e, long int n) {
    int encryptedNum = exponentiationRemainder(num, e, n);

    return encryptedNum;
}

__device__ int decrypt(long int num, long int e, long int phiN, long int n) {
    long int d;
    long int d1;
    long int gcd;
    gcd = extendedEuclidGcd(e, phiN, &d, &d1);
    int decryptedNum;
    if(d > 0) {
        decryptedNum = exponentiationRemainder(num, d, n);
    } else {
        // special case, use d + Phi-n as exponent
        decryptedNum = exponentiationRemainder(num, d + phiN, n);
    }

    return decryptedNum;
}

__global__ void run(int *deviceSieve, int compressedSize, int *deviceDecryption) {
    int idx = threadIdx.x + blockIdx.x*blockDim.x;

    // initialize state...used for curand
    curandState state;
    curand_init(clock64(), idx, 0, &state);

    int wordNum = generateWordNum(4, idx, state);

    // max and min index limits for sieve of eratosthenes
    int maxLimit = floorf(compressedSize/10);
    int minLimit = maxLimit - floorf(maxLimit/10);
    int randomNumberOne = getRandomNumber(idx, minLimit, maxLimit, state);
    int randomNumberTwo = getRandomNumber(idx, minLimit, maxLimit, state);
    while(randomNumberOne == randomNumberTwo) {
        randomNumberTwo = getRandomNumber(idx, minLimit, maxLimit, state);
    }
    int p = deviceSieve[randomNumberOne];
    int q = deviceSieve[randomNumberTwo];

    long int n = p*q;
    long int phiN = (p-1)*(q-1);
    long int e = generateE(idx, state, phiN);

    int encrypted = encrypt(wordNum, e, n);
    int decrypted = decrypt(encrypted, e, phiN, n);

    // keep the decryption idx value at 0 if encryption/decryption failed
    if(wordNum == decrypted) {
        deviceDecryption[idx] = wordNum;
    }
}

int main(void) {
    // get number of blocks for execution
    int numBlocks;
    cout << "Enter number of blocks: ";
    cin >> numBlocks;

    int maxPrime = 50000;
    int sieve[maxPrime];
    generateSieve(sieve, maxPrime);
    int compressedSize = getCompressedSize(sieve, maxPrime);
    int sieveOfEratosthenes[compressedSize];
    generateCompressedSieve(sieve, sieveOfEratosthenes, maxPrime);

    int *deviceSieve;
    cudaMalloc(&deviceSieve, sizeof(int)*compressedSize);
    cudaMemcpy(deviceSieve, sieveOfEratosthenes, sizeof(int)*compressedSize, cudaMemcpyHostToDevice);

    // create decryption array and initialize values to 0
    int *decryption = (int*)malloc(1024*numBlocks*sizeof(int));
    for(int i = 0; i < 1024*numBlocks; ++i) {
        decryption[i] = 0;
    }
    int *deviceDecryption;
    cudaMalloc(&deviceDecryption, sizeof(int)*1024*numBlocks);
    cudaMemcpy(deviceDecryption, decryption, sizeof(int)*1024*numBlocks, cudaMemcpyHostToDevice);

    auto start = chrono::high_resolution_clock::now();
    run<<<numBlocks,1024>>>(deviceSieve, compressedSize, deviceDecryption);
    cudaDeviceSynchronize();
    auto end = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::milliseconds>(end-start);

    cudaMemcpy(decryption, deviceDecryption, sizeof(int)*1024*numBlocks, cudaMemcpyDeviceToHost);

    // count successes and failures
    int success = 0;
    int failure = 0;
    for(int j = 0; j < 1024*numBlocks; ++j) {
        if(decryption[j]) {
            success++;
        } else {
            failure++;
        }
    }

    // display results to user
    cout << "SUCCESS: " + to_string(success) + "\n";
    cout << "FAILURE: " + to_string(failure) + "\n";
    cout << "TOTAL TIME: " + to_string((float)duration.count()) + " MILLISECONDS" + "\n";

    // free memory at end of execution
    cudaFree(deviceSieve);
    cudaFree(deviceDecryption);
    free(decryption);

    return 1;
}