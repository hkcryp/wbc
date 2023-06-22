#include <iostream>
#include <cstdint>
#include <fstream>
#include <cstring>
#include <ctime>
#include <assert.h>
#include <x86intrin.h>

#define blocknum 96
#define readSize  0x400000000

#define MAX_BRANCHES 8
#define SPARKLE_STATE   256
#define STATE_BRANS (SPARKLE_STATE/64)
#define STATE_WORDS (SPARKLE_STATE/32)
#define round 14

#define ROT(x, n) (((x) >> (n)) | ((x) << (32-(n))))
#define ELL(x) (ROT(((x) ^ ((x) << 16)), 16))

#pragma intrinsic(__rdtsc)
uint64_t start_rdtsc(){return __rdtsc();}
uint64_t end_rdtsc(){return __rdtsc();}

using namespace std;

typedef struct {
  uint32_t x[MAX_BRANCHES];
  uint32_t y[MAX_BRANCHES];
} SparkleState;

uint32_t RCON[MAX_BRANCHES] = {0xB7E15162, 0xBF715880, 0x38B4DA56, 0x324E7738, 0xBB1185EB, 0x4F7C7B57, 0xCFBFA1C8, 0xC2B3293D};

void nonlinear_layer(SparkleState& state, char* memblock){
  uint64_t tmp;

  for(int j=0;j<STATE_BRANS;j++){
    tmp=(uint64_t)4*state.x[j];
    memcpy(&state.x[j],&memblock[tmp],sizeof(int));
    tmp=(uint64_t)4*state.y[j];
    memcpy(&state.y[j],&memblock[tmp],sizeof(int));
  }
}

void linear_layer(SparkleState& state)
{
  int i, b = STATE_BRANS/2;
  uint32_t *x = state.x, *y = state.y;
  uint32_t tmp;
  
  tmp = 0;
  for(i = 0; i < b; i++){tmp ^= x[i];}
  tmp = ELL(tmp);
  for(i = 0; i < b; i ++){y[i+b] ^= (tmp ^ y[i]);}
  
  tmp = 0;
  for(i = 0; i < b; i++){tmp ^= y[i];}
  tmp = ELL(tmp);
  for(i = 0; i < b; i ++){x[i+b] ^= (tmp ^ x[i]);}
  
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

void sparkle_ref(SparkleState& state, char* memblock)
{
  int i, j;

    uint64_t tmp;  
    for(i = 0; i < round; i++) {
      state.y[0] ^= RCON[i%MAX_BRANCHES];
      state.y[1] ^= i;
      nonlinear_layer(state, memblock);
      linear_layer(state);
    }
}

void generatemessage(uint32_t input[blocknum][STATE_WORDS]){
	for(int i=0; i<blocknum; i++){
    for(int j=0;j<STATE_WORDS;j++){
      input[i][j]=rand();
    }
	}
}

int main(){
  uint32_t input[blocknum][STATE_WORDS]={};
  SparkleState state[blocknum] = {};

  int iter=100000;
  
  uint64_t start_time, end_time, cyc=0;
  srand(time(0));

  ifstream table("~/table/sbox32.bin", ios::in | ios::binary );

  if(table.is_open()){
    char* memblock = new char[readSize];
    table.read(memblock,readSize);
    table.close(); 

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
        sparkle_ref(state[j], memblock);
        end_time=_rdtsc();
        cyc+=end_time-start_time;
      }
    }

    delete[] memblock;
  }

  uint64_t aver=cyc/iter;

  cout << "average cost " << aver << " CPU cycles for " << iter << " test" << endl;
	uint64_t cpb = aver/(4*STATE_WORDS*blocknum);
  cout << "average CPB for decryption is " << cpb << endl;

  return 0;
} 