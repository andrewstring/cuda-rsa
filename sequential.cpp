#include <iostream>
#include <chrono>
#include <cmath>

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

// compressed size will be used in main function
int getCompressedSize(int *sieve, int maxPrime) {
    int compressedSize = 0;
    for(int i = 0; i < maxPrime; ++i) {
        if(sieve[i] != 0) {
            compressedSize++;
        }
    }

    return compressedSize - 1;
}

int getRandomNumber(int min, int max) {
    int difference = max - min + 1;
    return rand() % (difference) + min;
}

// generate numerical representation of word
int generateWordNum(int size) {
    int wordNum = 0;
    int counter = 0;
    for(int j = 0; j < size; ++j) {
        int randomNum = getRandomNumber(1, 26);
        wordNum = wordNum + randomNum * pow(27, counter);
        counter++;
    }

    return wordNum;
}

int extendedEuclidGcd(long int first, long int second, long int *x, long int *y) {
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

long int euclidGcd(long int first, long int second) {
    if(!second) {
        return first;
    } else {
        long int remainder = first%second;
        return euclidGcd(second, remainder);
    }
}

// e is an odd number that is relatively prime to Phi-n
int generateE(long int phiN) {
    int e = 2;
    while(e == 2 || e%2 == 0 || euclidGcd(e, phiN) != 1) {
        e = getRandomNumber(100, 500);
    }

    return e;
}

// calculates modulus of exponentiation without overflow
int exponentiationRemainder(long int a, long int b, long int c) {
    int result = 1;
    while(b > 0) {
        if(b%2 == 1) {
            result = (result*a) % c;
        }
        b = floorf(b/2);
        a = (a*a)%c;
    }

    return result;
}

int encrypt(long int num, long int e, long int n) {
    int encryptedNum = exponentiationRemainder(num, e, n);

    return encryptedNum;
}

int decrypt(long int num, long int e, long int phiN, long int n) {
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

void run(int idx, int *sieve, int compressedSize, int *decryption) {
    // use time as seed
    srand(time(NULL));

    int wordNum = generateWordNum(4);

    // max and min index limits for sieve of eratosthenes
    int maxLimit = floor(compressedSize/10);
    int minLimit = maxLimit - floor(maxLimit/10);
    int p = sieve[getRandomNumber(minLimit, maxLimit)];
    int q = sieve[getRandomNumber(minLimit, maxLimit)];
    while(p == q) {
        q = sieve[getRandomNumber(minLimit, maxLimit)];
    }

    long int n = p*q;
    long int phiN = (p-1)*(q-1);
    long int e = generateE(phiN);

    int encrypted = encrypt(wordNum, e, n);
    int decrypted = decrypt(encrypted, e, phiN, n);

    // keep the decryption idx value at 0 if encryption/decryption failed
    if(wordNum == decrypted) {
        decryption[idx] = wordNum;
    }
}

int main(void) {
    // get number multiple of 1024 for program execution
    int numMultiple;
    cout << "Enter multiple of 1024 encryptions/decryptions to run: ";
    cin >> numMultiple;

    int maxPrime = 50000;
    int sieve[maxPrime];
    generateSieve(sieve, maxPrime);
    int compressedSize = getCompressedSize(sieve, maxPrime);
    int sieveOfEratosthenes[compressedSize];
    generateCompressedSieve(sieve, sieveOfEratosthenes, maxPrime);

    // create decryption array and initialize values to 0
    int decryption[numMultiple*1024];
    for(int i = 0; i < 1024*numMultiple; ++i) {
        decryption[i] = 0;
    }

    auto start = chrono::high_resolution_clock::now();
    for(int j = 0; j < 1024*numMultiple; ++j) {
        run(j, sieveOfEratosthenes, compressedSize, decryption);
    }
    auto end = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::milliseconds>(end-start);

    // count successes and failures
    int success = 0;
    int failure = 0;
    for(int k = 0; k < 1024*numMultiple; ++k) {
        if(decryption[k]) {
            success++;
        } else {
            failure++;
        }
    }

    // display results to user
    cout << "SUCCESS: " + to_string(success) + "\n";
    cout << "FAILURE: " + to_string(failure) + "\n";
    cout << "TOTAL TIME: " + to_string((float)duration.count()) + " MILLISECONDS" + "\n";

    return 1;
}