#include <iostream>
#include <cstdint>
#include <fstream>
#include <cstring>
#include <ctime>
#include <assert.h>
#include <x86intrin.h>
#include "tbox.h"

#define blocknum 96
#define readSize  0x400000000

#define MAX_BRANCHES 8 // 8
#define SPARKLE_STATE   256
#define STATE_BRANS (SPARKLE_STATE/32)
#define STATE_WORDS (SPARKLE_STATE/16)
#define round 8

#define ROT(x, n) (((x) >> (n)) | ((x) << (16-(n))))
#define ELL(x) (ROT(((x) ^ ((x) << 8)), 8))

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

void linear_layer(SparkleState& state, int brans)
{
  int i, b = brans/2;
  uint16_t *x = state.x, *y = state.y;
  uint16_t tmp=0;
  
  tmp = 0;
  for(i = 0; i < b; i++){tmp ^= x[i];}
    
  tmp = ELL(tmp); 
  for(i = 0; i < b; i ++){y[i+b] ^= (tmp ^ y[i]); }
  
  tmp = 0;
  for(i = 0; i < b; i++){tmp ^= y[i]; }
    
  tmp = ELL(tmp); 
  for(i = 0; i < b; i ++){x[i+b] ^= (tmp ^ x[i]); }
    
  tmp = x[0];
  for (i = 0; i < b - 1; i++) {
    x[i] = x[i+b+1];
    x[i+b+1] = x[i+1];
  }
  x[b-1] = x[b];
  x[b] = tmp;
  
  tmp = y[0];
  for (i = 0; i < b - 1; i++) {
    y[i] = y[i+b+1];
    y[i+b+1] = y[i+1];
  }
  y[b-1] = y[b];
  y[b] = tmp;
}

void sparkle_ref(SparkleState& state, int brans)
{
  int i, j;

  assert(((brans & 1) == 0) && (brans >= 4) && (brans <= MAX_BRANCHES));

  uint32_t tmp;  
  for(i = 0; i < round; i++) {
    state.y[0] ^= RCON[i%MAX_BRANCHES]; 
    state.y[1] ^= i; 
    for(j=0;j<brans;j++){
      state.x[j] = invT[state.x[j]];
      state.y[j] = invT[state.y[j]];
    }
    linear_layer(state, brans);
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
      start_time =_rdtsc();        
      sparkle_ref(state[j], STATE_BRANS);
      end_time =_rdtsc();
      cyc += end_time - start_time;
    }

  }

  uint64_t aver = cyc/iter;

  cout << "average cost " << aver << " CPU cycles for " << iter << " test" << endl;
	uint64_t cpb = aver/(2*STATE_WORDS*blocknum);
  cout << "average CPB for decryption is " << cpb << endl;

  return 0;
} 
