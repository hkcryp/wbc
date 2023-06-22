#include <iostream>
#include <cstdint>
#include <fstream>
#include <cstring>
#include <ctime>
#include <assert.h>
#include <x86intrin.h>

#define blocknum 96

#define round1632 24
#define ROTL8(x,r)  (((x)<<(r)) | (x>>(8-(r))))
#define ROTR8(x,r) (((x)>>(r)) | ((x)<<(8-(r))))
#define newER8(x,y,k_left,k_right) (x^=k_left, y^=k_right, x=ROTR8(x,7), x+=y, y=ROTL8(y,2), y^=x)
#define newDR8(x,y,k_left,k_right) (y^=x, y=ROTR8(y,2), x-=y, x=ROTL8(x,7), y^=k_right, x^=k_left)
#define ROT(x, n) (((x) >> (n)) | ((x) << (16-(n))))
#define ELL(x) (ROT(((x) ^ ((x) << 8)), 8))

#define MAX_BRANCHES 8
#define SPARKLE_STATE   256
#define STATE_BRANS (SPARKLE_STATE/32)
#define STATE_WORDS (SPARKLE_STATE/16)

#define round 8

#pragma intrinsic(__rdtsc)
uint64_t start_rdtsc()
{
	return __rdtsc();
}
uint64_t end_rdtsc()
{
	return __rdtsc();
}

using namespace std;

typedef struct {
  uint16_t x[MAX_BRANCHES];
  uint16_t y[MAX_BRANCHES];
} SparkleState;

static const uint16_t RCON[MAX_BRANCHES] = { 0xB7E1, 0x5162, 0xBF71, 0x5880, 0x38B4, 0xDA56, 0x324E, 0x7738 };

void newSpeck1632Encrypt(uint8_t Pt[], uint8_t Ct[], uint8_t rk[])
{
	uint8_t i;
	Ct[0] = Pt[0]; Ct[1] = Pt[1];
	for (i = 0; i < 2 * round1632; i += 2)
		newER8(Ct[1], Ct[0], rk[i], rk[i + 1]);
}

void lookup_layer(SparkleState& state, int brans, uint8_t roundkey[]){

  uint8_t arrPx[brans][2]={}, arrCx[brans][2]={};

  for(int i=0; i<brans; i++){
    arrPx[i][0] = (uint8_t)(state.x[i] >> 8);
    arrPx[i][1] = (uint8_t)(state.x[i]);
    newSpeck1632Encrypt(arrPx[i], arrCx[i], roundkey);
    state.x[i] = (arrCx[i][0] << 8) ^ arrCx[i][1];    
  }

  uint8_t arrPy[brans][2]={}, arrCy[brans][2]={};
  
  for(int i=0; i<brans; i++){
    arrPy[i][0] = (uint8_t)(state.y[i] >> 8);
    arrPy[i][1] = (uint8_t)(state.y[i]);
    newSpeck1632Encrypt(arrPy[i], arrCy[i], roundkey);
    state.y[i] = (arrCy[i][0] << 8) ^ arrCy[i][1];    
  }
}

void linear_layer_inv(SparkleState& state, int brans)
{
  int i, b = brans/2;
  uint16_t *x = state.x, *y = state.y;
  uint16_t tmp=0;
  
  tmp = x[b-1];
  for (i = b - 1; i > 0; i--) {
    x[i] = x[i+b];
    x[i+b] = x[i-1];
  }
  x[0] = x[b];
  x[b] = tmp;
  
  tmp = y[b-1];
  for (i = b - 1; i > 0; i--) {
    y[i] = y[i+b];
    y[i+b] = y[i-1];
  }
  y[0] = y[b];
  y[b] = tmp;
  
  tmp = 0;
  for(i = 0; i < b; i ++)
    tmp ^= y[i];
  tmp = ELL(tmp);
  for(i = 0; i < b; i ++)
    x[i+b] ^= (tmp ^ x[i]);
  
  tmp = 0;
  for(i = 0; i < b; i ++)
    tmp ^= x[i];
  tmp = ELL(tmp);
  for(i = 0; i < b; i ++)
    y[i+b] ^= (tmp ^ y[i]);
}


void sparkle_inv_ref(SparkleState& state, int brans)
{
  int i, j;  // Step and branch counter
  uint8_t rkeys[48]={0x46,0xb9,0xdd,0x2b,0x0b,0xa8,0x8d,0x13,0x23,0x3b,0x3f,0xeb,0x74,0x3e,0xeb,0x24,0x3f,0xcd,0x52,0xea,0x62,0xb8,0x1b,0x82,0xb5,0x0c,0x27,0x64,0x6e,0xd5,0x76,0x2f,0xd7,0x5d,0xc4,0xdd,0xd8,0xc0,0xf2,0x00,0xcb,0x05,0x01,0x9d,0x67,0xb5,0x92,0xf6};
  
  assert(((brans & 1) == 0) && (brans >= 4) && (brans <= MAX_BRANCHES));
  
  for(i = round - 1; i >= 0; i--) {
    linear_layer_inv(state, brans);
    lookup_layer(state, brans, rkeys);
    state.y[1] ^= i;
    state.y[0] ^= RCON[i%MAX_BRANCHES];
  }
}

void generatemessage(uint16_t input[blocknum][STATE_WORDS]){
	for(int i=0; i<blocknum; i++){
    for(int j=0;j<STATE_WORDS;j++){
      input[i][j]=rand();
    }
	}
}

int main(){
  uint16_t input[blocknum][STATE_WORDS]={};
  SparkleState state[blocknum] = {};

  int iter=100000;
  
  uint64_t start_time, end_time, cyc=0;
  srand(time(0));

  generatemessage(input); 

  for(int count=0; count<iter; count ++){
    
    for(int j=0; j<blocknum; j++){      
      for (int i = 0; i < (STATE_WORDS)/2; i++) {
        state[j].x[i] = input[j][2*i];
        state[j].y[i] = input[j][2*i+1];
      }
    }
    
    for(int j=0; j<blocknum; j++){
      start_time=_rdtsc();
      sparkle_inv_ref(state[j], STATE_BRANS);
      end_time=_rdtsc();
      cyc+=end_time-start_time;
    }

  }

  uint64_t aver=cyc/iter;

  cout << "average cost " << aver << " CPU cycles for " << iter << " test" << endl;
	uint64_t cpb = aver/(2*STATE_WORDS*blocknum);
  cout << "average CPB for decryption is " << cpb << endl;

  return 0;
} 

