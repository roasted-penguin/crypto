/*
 * Copyright 2020-2022. Heekuck Oh, all rights reserved
 * 이 프로그램은 한양대학교 ERICA 소프트웨어학부 재학생을 위한 교육용으로 제작되었다.
 */
#ifdef __linux__
#include <bsd/stdlib.h>
#elif __APPLE__
#include <stdlib.h>
#else
#include <stdlib.h>
#endif
#include "mRSA.h"

/*
 * mod_add() - computes a + b mod m
 */
static uint64_t mod_add(uint64_t a, uint64_t b, uint64_t m)
{
  if(a>m) a = a%m;
  if(b>m) b = b%m;
  if(a>=m-b){
    return a-(m-b);
  }
  else return a+b;
}

/*
 * mod_mul() - computes a * b mod m
 */
static uint64_t mod_mul(uint64_t a, uint64_t b, uint64_t m)
{
  uint64_t r = 0;
  while(b>0){
    if(b&1){
      r = mod_add(r,a,m);//r+=a
    }
    b = b >> 1;//b/=2
    a = mod_add(a,a,m);//a*=2
  }
  return r;
}

/*
 * mod_pow() - computes a^b mod m
 */
static uint64_t mod_pow(uint64_t a, uint64_t b, uint64_t m)
{
  a = a%m;
  uint64_t r = 1;
  while(b>0){
    if (b&1){
      r = mod_mul(r,a,m);
    }
    b = b>>1;
    a = mod_mul(a,a,m);
  }
  return r;
}

/*
 * gcd() - Euclidean algorithm
 */
static uint64_t gcd(uint64_t a, uint64_t b)
{
  if(a==0){//a==0일때 선처리
    return b;
  }
  else{
    uint64_t n;
    while(b!=0){
      n = a%b;
      a = b;
      b = n;
    }
    return a;
  }
}

static uint64_t carmichael(uint64_t a, uint64_t b){
  uint64_t n = a-1;
  n = n/gcd(a-1,b-1);
  n = n * (b-1);
  return n;
}
/*
 * mul_inv() - computes multiplicative inverse a^-1 mod m
 * It returns 0 if no inverse exist.
 */
 static uint64_t umul_inv(uint64_t a, uint64_t m)
 {
   uint64_t d0, d1, tmp1;//d는 항상 양의 정수
   int64_t x0, x1, tmp2;//x는 정수, 즉 음의 정수가 될 수 있음
   d0 = a; d1 = m; x0 = 1; x1 = 0;
   uint64_t q;
   while(d1 > 1){
     q = d0/d1;
     tmp1 = d0 - q*d1;//d1이 될 수 연산
     d0 = d1;//swap
     d1 = tmp1;

     tmp2 = x0 - q*x1;//x1이 될 수 연산
     x0 = x1;//swap
     x1 = tmp2;
   }
   if(d1 == 1) return (x1 > 0 ? x1 : x1+m);
   else return 0;
 }

/*
 * Miller-Rabin Primality Testing against small sets of bases
 *
 * if n < 2^64,
 * it is enough to test a = 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, and 37.
 *
 * if n < 3317044064679887385961981,
 * it is enough to test a = 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, and 41.
 */
static const uint64_t a_arr[BASELEN] = {2,3,5,7,11,13,17,19,23,29,31,37};

/*
 * miller_rabin() - Miller-Rabin Primality Test (deterministic version)
 *
 * n > 3, an odd integer to be tested for primality
 * It returns 1 if n is prime, 0 otherwise.
 */
static int miller_rabin(uint64_t n)
{
  if(n==2) return PRIME;
  if(n==1 || n%2==0) return COMPOSITE;
  else{
    //printf("%d in 1\n",n);
    //n-1 = 2^k*q형태로 변환.
    uint64_t q = n-1;
    int k = 0;
    while(q%2==0){
      q = q>>1;
      k++;
    }

    int cnt = 0;
    for(int i=0;i<12;i++){//a크기=12만큼 반복
      uint64_t a = a_arr[i];
      if(n==a){
        return PRIME;
      }
      uint64_t x = mod_pow(a,q,n);
      if(x==1){
        cnt++; continue;
      }
      if(cnt<=i){
        for(int j=0;j<=k-1;j++){
          if(x==n-1) {
            cnt++; break;
          }
          x = mod_mul(x,x,n);
        }
      }
      if(cnt<=i) return COMPOSITE;//cnt안걸렸으면 composite 리턴
    }
    return PRIME;
  }
}

/*
 * mRSA_generate_key() - generates mini RSA keys e, d and n
 *
 * Carmichael's totient function Lambda(n) is used.
 */
void mRSA_generate_key(uint64_t *e, uint64_t *d, uint64_t *n)
{
  uint64_t p,q;
  //소수 2개 p,q 생성 : p와 q의 길이가 비슷할수록 더 안전하다는 점을 참고.
  do{
    p = arc4random_uniform(0x80000000)+0x80000000;
    if(p%2==0) p = p+1;
  } while(!miller_rabin(p));

  uint64_t qmin = (0x8000000000000000/p) + 1;
  uint64_t qmax = 0xffffffffffffffff/p;
  do{
    q = arc4random_uniform(qmax-qmin+1)+qmin;
    if(q%2==0) q = q+1;
  } while(!miller_rabin(q));
  *n = p*q;

  uint64_t cmn = carmichael(p,q);//'C'ar'M'ichael 'N'

  //*e = 65537; //make encryption faster
  do{
    *e = arc4random_uniform(cmn-4)+4;//3<e<𝜆(n) 3은 취약하므로 피한다
  } while(gcd(*e,cmn)!=1);

  *d = umul_inv(*e,cmn);
}

/*
 * mRSA_cipher() - compute m^k mod n
 *
 * If data >= n then returns 1 (error), otherwise 0 (success).
 */
int mRSA_cipher(uint64_t *m, uint64_t k, uint64_t n)
{
  if(*m>=n) return 1;
  *m = mod_pow(*m,k,n);
  return 0;
}
