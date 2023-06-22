#include <iostream>
#include <cstdint>
#include <fstream>
#include <cstring>
#include <ctime>
#include <assert.h>
#include <x86intrin.h>

#define blocknum 96

#define MAX_BRANCHES 8
#define SPARKLE_STATE   256
#define STATE_BRANS (SPARKLE_STATE/64)
#define STATE_WORDS (SPARKLE_STATE/32)
#define round 14

#define ROT(x, n) (((x) >> (n)) | ((x) << (32-(n))))
#define ELL(x) (ROT(((x) ^ ((x) << 16)), 16))
#define xtime(a) (((a) << (1)) ^ ((((a) >> (7)) & 1) * (0x1b)))

#pragma intrinsic(__rdtsc)
uint64_t start_rdtsc(){return __rdtsc();}
uint64_t end_rdtsc(){return __rdtsc();}

using namespace std;

typedef struct {
  uint32_t x[MAX_BRANCHES];
  uint32_t y[MAX_BRANCHES];
} SparkleState;

uint32_t RCON[MAX_BRANCHES] = { 0xB7E15162, 0xBF715880, 0x38B4DA56, 0x324E7738, 0xBB1185EB, 0x4F7C7B57, 0xCFBFA1C8, 0xC2B3293D };

uint8_t invSB[256]={0x52, 0x09, 0x6a, 0xd5, 0x30, 0x36, 0xa5, 0x38, 0xbf, 0x40, 0xa3, 0x9e, 0x81, 0xf3, 0xd7, 0xfb,
		            0x7c, 0xe3, 0x39, 0x82, 0x9b, 0x2f, 0xff, 0x87, 0x34, 0x8e, 0x43, 0x44, 0xc4, 0xde, 0xe9, 0xcb,
		            0x54, 0x7b, 0x94, 0x32, 0xa6, 0xc2, 0x23, 0x3d, 0xee, 0x4c, 0x95, 0x0b, 0x42, 0xfa, 0xc3, 0x4e,
		            0x08, 0x2e, 0xa1, 0x66, 0x28, 0xd9, 0x24, 0xb2, 0x76, 0x5b, 0xa2, 0x49, 0x6d, 0x8b, 0xd1, 0x25,
		            0x72, 0xf8, 0xf6, 0x64, 0x86, 0x68, 0x98, 0x16, 0xd4, 0xa4, 0x5c, 0xcc, 0x5d, 0x65, 0xb6, 0x92,
		            0x6c, 0x70, 0x48, 0x50, 0xfd, 0xed, 0xb9, 0xda, 0x5e, 0x15, 0x46, 0x57, 0xa7, 0x8d, 0x9d, 0x84,
		            0x90, 0xd8, 0xab, 0x00, 0x8c, 0xbc, 0xd3, 0x0a, 0xf7, 0xe4, 0x58, 0x05, 0xb8, 0xb3, 0x45, 0x06,
		            0xd0, 0x2c, 0x1e, 0x8f, 0xca, 0x3f, 0x0f, 0x02, 0xc1, 0xaf, 0xbd, 0x03, 0x01, 0x13, 0x8a, 0x6b,
		            0x3a, 0x91, 0x11, 0x41, 0x4f, 0x67, 0xdc, 0xea, 0x97, 0xf2, 0xcf, 0xce, 0xf0, 0xb4, 0xe6, 0x73,
		            0x96, 0xac, 0x74, 0x22, 0xe7, 0xad, 0x35, 0x85, 0xe2, 0xf9, 0x37, 0xe8, 0x1c, 0x75, 0xdf, 0x6e,
		            0x47, 0xf1, 0x1a, 0x71, 0x1d, 0x29, 0xc5, 0x89, 0x6f, 0xb7, 0x62, 0x0e, 0xaa, 0x18, 0xbe, 0x1b,
		            0xfc, 0x56, 0x3e, 0x4b, 0xc6, 0xd2, 0x79, 0x20, 0x9a, 0xdb, 0xc0, 0xfe, 0x78, 0xcd, 0x5a, 0xf4,
		            0x1f, 0xdd, 0xa8, 0x33, 0x88, 0x07, 0xc7, 0x31, 0xb1, 0x12, 0x10, 0x59, 0x27, 0x80, 0xec, 0x5f,
		            0x60, 0x51, 0x7f, 0xa9, 0x19, 0xb5, 0x4a, 0x0d, 0x2d, 0xe5, 0x7a, 0x9f, 0x93, 0xc9, 0x9c, 0xef,
		            0xa0, 0xe0, 0x3b, 0x4d, 0xae, 0x2a, 0xf5, 0xb0, 0xc8, 0xeb, 0xbb, 0x3c, 0x83, 0x53, 0x99, 0x61,
		            0x17, 0x2b, 0x04, 0x7e, 0xba, 0x77, 0xd6, 0x26, 0xe1, 0x69, 0x14, 0x63, 0x55, 0x21, 0x0c, 0x7d};

void invSboxLayer(uint32_t& state){
  uint32_t tmp=0;

  for(int j=3; j>=0; j--){tmp|=(invSB[(state>>8*j)&0xFF]<<8*j);}

  state=tmp; 
}

void addroundkeyLayer(uint32_t& state, uint32_t roundkey){state^=roundkey;}

void invmixColumnLayer(uint32_t &state){

    uint8_t tmp[4] = {};
    uint8_t st[4] = {};

    for(int j=0;j<4;j++){
        tmp[j] = state&0xFF; 
        state=state>>8;
    }
    
    st[0] = xtime(xtime(xtime(tmp[0]^tmp[1]^tmp[2]^tmp[3]))) ^ xtime(xtime(tmp[0]^tmp[2])) ^ xtime(tmp[0]^tmp[1]) ^ tmp[1]^tmp[2]^tmp[3];
    st[1] = xtime(xtime(xtime(tmp[0]^tmp[1]^tmp[2]^tmp[3]))) ^ xtime(xtime(tmp[1]^tmp[3])) ^ xtime(tmp[1]^tmp[2]) ^ tmp[0]^tmp[2]^tmp[3];
    st[2] = xtime(xtime(xtime(tmp[0]^tmp[1]^tmp[2]^tmp[3]))) ^ xtime(xtime(tmp[0]^tmp[2])) ^ xtime(tmp[2]^tmp[3]) ^ tmp[0]^tmp[1]^tmp[3];
    st[3] = xtime(xtime(xtime(tmp[0]^tmp[1]^tmp[2]^tmp[3]))) ^ xtime(xtime(tmp[1]^tmp[3])) ^ xtime(tmp[0]^tmp[3]) ^ tmp[0]^tmp[1]^tmp[2];
    
    for(int j=3;j>=0;--j){
        state = state<<8; 
        state ^= st[j];
    } 
}

void lookup_layer(SparkleState &state){

  uint32_t key[0x10]={0x46b9dd2b,0x0ba88d13,0x233b3feb,0x743eeb24,0x3fcd52ea,0x62b81b82,0xb50c2764,0x6ed5762f,
                      0xd75dc4dd,0xd8c0f200,0xcb05019d,0x67b592f6,0xfc821c49,0x479ab486,0x40292eac,0xb3b7c4be};

  for(int j=0; j<STATE_BRANS; j++){
    for(int i=15; i>0; i--){
      addroundkeyLayer(state.x[j],key[i]);
      invmixColumnLayer(state.x[j]);
      invSboxLayer(state.x[j]);
    }
    addroundkeyLayer(state.x[j], key[0]);

    for(int i=15; i>0; i--){
      addroundkeyLayer(state.y[j],key[i]);
      invmixColumnLayer(state.y[j]);
      invSboxLayer(state.y[j]);
    }
    addroundkeyLayer(state.y[j], key[0]);
  }

}

void linear_layer_inv(SparkleState& state)
{
  int i, b = STATE_BRANS/2;
  uint32_t *x = state.x, *y = state.y;
  uint32_t tmp;
  
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
  for(i = 0; i < b; i ++){tmp ^= y[i];}    
  tmp = ELL(tmp);
  for(i = 0; i < b; i ++){x[i+b] ^= (tmp ^ x[i]);}
    
  tmp = 0;
  for(i = 0; i < b; i ++){tmp ^= x[i];} 
  tmp = ELL(tmp);
  for(i = 0; i < b; i ++){y[i+b] ^= (tmp ^ y[i]);}
}

void sparkle_inv_ref(SparkleState& state)
{  
  for(int i = round - 1; i >= 0; i--) {
    linear_layer_inv(state);
    lookup_layer(state);
    state.y[1] ^= i;
    state.y[0] ^= RCON[i%MAX_BRANCHES];
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
      sparkle_inv_ref(state[j]);
      end_time=_rdtsc();
      cyc+=end_time-start_time;
    }
  }

  uint64_t aver=cyc/iter;

  cout << "average cost " << aver << " CPU cycles for " << iter << " test" << endl;
	uint64_t cpb = aver/(4*STATE_WORDS*blocknum);
  cout << "average CPB for decryption is " << cpb << endl;

  return 0;
} 
