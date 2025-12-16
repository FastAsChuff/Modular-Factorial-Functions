//===============================================================================================================================
// kfactorialmodn.c: START
//===============================================================================================================================
// Author: Simon Goater Dec 2025
// -----------------------------
// Simple functions to calculate k! mod n.
// getkfactmodn1() - k,n < 2^64 (Slowest on x64)
// Define macro KFACTORIALMODUSEMG to use the following functions.
//    Require https://github.com/FastAsChuff/Fast-Modular-Exponentiation/blob/main/modpowu64.c
//    and https://github.com/FastAsChuff/Fast-Modular-Inverse-Modulo-Powers-Of-2/blob/main/fastmodinvpow2fns.c
// getkfactmodn6() - n < 2^64, k < 2^32, n must be odd and pre-computed array primes must contain all primes <= k in ascending 
//    order. (Fastest on x64)
// getkfactmodn3() - k,n < 2^64, n must be odd. 
// COPYRIGHT: Software is given as is without guarantee or warranty. It is offered free to use, modify, copy, or distribute
// with conspicuous attribution for any purpose.
//===============================================================================================================================

#ifndef U128
  #define U128 unsigned __int128
#endif

#ifdef KFACTORIALMODUSEMG

uint64_t getkfactmodn3(uint64_t k, uint64_t n) {
  if (k >= n) return 0;
  uint64_t twoto64modn = 0;
  uint64_t ninv = modinv64x(n);
  uint64_t r = tou64mg(1, n, &twoto64modn);
  uint64_t fr1 = r;
  uint64_t ir1 = r;
  uint64_t kk = k/2;
  uint64_t ir2 = tou64mg(kk+1, n, &twoto64modn);
  uint64_t fr2 = r;
  for (uint64_t i=1; i<=kk; i++) {
    fr1 = modprodu64mg(fr1, ir1, n, ninv, twoto64modn);
    fr2 = modprodu64mg(fr2, ir2, n, ninv, twoto64modn);
    ir1 = modsumu64mg(ir1, r, n, ninv, twoto64modn);
    ir2 = modsumu64mg(ir2, r, n, ninv, twoto64modn);
  }
  uint64_t f1 = fromu64mg(fr1, n, ninv, twoto64modn);
  uint64_t f2 = fromu64mg(fr2, n, ninv, twoto64modn);
  if (k & 0x1u) f1 = ((U128)k*f1) % n;
  return ((U128)f2*f1) % n;
}

uint64_t getkfactmodn6(uint32_t k, uint64_t n, uint32_t numprimes, uint32_t *primes) {
  // Assumes primes contains all primes <= k.
  // assert(n & 0x1u);
  if (k >= n) return 0;
  uint32_t pix = 0;
  uint64_t res = 1;
  uint64_t ei = 0;
  uint64_t eiprev = 0;
  uint64_t p = 2;
  uint64_t prodp = 1; // res *= prodp^eiprev mod n
  while (true) { // Process primes with common ei
    while (true) {
      uint64_t primepower = primes[pix];
      ei = 0;
      while (true) { // Calculate ei
        uint64_t term = (k/primepower);
        ei += term;
        primepower = (primepower*primes[pix]);
        if (primepower > k) break;
      }
      if (eiprev == 0) eiprev = ei;
      if (ei != eiprev) break;
      prodp = ((U128)prodp * p) % n;
      pix++;
      if (pix >= numprimes) break;
      if (primes[pix] > k) break;
      p = primes[pix];
    }
    prodp = modpowu64(prodp, eiprev, n);
    res = ((U128)res * prodp) % n;
    if (pix >= numprimes) break;
    if (primes[pix] > k) break;
    prodp = 1;
    eiprev = ei;
  }
  return res;
}
#endif

uint64_t getkfactmodn1(uint64_t k, uint64_t n) {
  if (k >= n) return 0;
  if (n >> 32) {
    uint64_t f = 1;
    for (uint64_t i=2; i<=k; i++) f = ((U128)f*i) % n;
    return f;
  } 
  uint32_t f = 1;
  for (uint32_t i=2; i<=k; i++) f = ((uint64_t)f*i) % n;
  return f;
}


//===============================================================================================================================
// kfactorialmodn.c: END
//===============================================================================================================================
