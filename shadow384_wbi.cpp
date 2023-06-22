#include <iostream>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <x86intrin.h>

#define blocknum 64
#define readSize  0x400000000

#define SHADOW_NS 15
#define SHADOW_NR 2 * SHADOW_NS
#define LS_ROWS 4
#define SHADOW_NBYTES 12
#define MLS_BUNDLES (SHADOW_NBYTES / LS_ROWS)

#define rotr(x,c) (((x) >> (c)) | ((x) << ((32) - (c))))

#pragma intrinsic(__rdtsc)
uint64_t start_rdtsc(){return __rdtsc();}
uint64_t end_rdtsc(){return __rdtsc();}

using namespace std;

static const uint32_t rc[SHADOW_NR][LS_ROWS] = {
  { 1, 0, 0, 0 }, // 0
  { 0, 1, 0, 0 }, // 1
  { 0, 0, 1, 0 }, // 2
  { 0, 0, 0, 1 }, // 3
  { 1, 1, 0, 0 }, // 4
  { 0, 1, 1, 0 }, // 5
  { 0, 0, 1, 1 }, // 6
  { 1, 1, 0, 1 }, // 7
  { 1, 0, 1, 0 }, // 8
  { 0, 1, 0, 1 }, // 9
  { 1, 1, 1, 0 }, // 10
  { 0, 1, 1, 1 }, // 11
  { 1, 0, 0, 0 }, // 12
  { 0, 1, 0, 0 }, // 13
  { 0, 0, 1, 0 }, // 14
  { 0, 0, 0, 1 }, // 15
  { 1, 1, 0, 0 }, // 16
  { 0, 1, 1, 0 }, // 17
  { 0, 0, 1, 1 }, // 18
  { 1, 1, 0, 1 }, // 19
  { 1, 0, 1, 0 }, // 20
  { 0, 1, 0, 1 }, // 21
  { 1, 1, 1, 0 }, // 22
  { 0, 1, 1, 1 }, // 23
  { 1, 0, 0, 0 }, // 24
  { 0, 1, 0, 0 }, // 25
  { 0, 0, 1, 0 }, // 26
  { 0, 0, 0, 1 }, // 27
  { 1, 1, 0, 0 }, // 28
  { 0, 1, 1, 0 }, // 29
};

void t_layer(uint32_t* state, char* memblock){
  uint64_t tmp;

  for(int i=0; i<LS_ROWS; i++){
    tmp=(uint64_t)4*state[i]; 
    memcpy(&state[i],&memblock[tmp],sizeof(int));
  }
}

void lbox(uint32_t* x, uint32_t* y) {
  uint32_t a, b, c, d;
  a = *x ^ rotr(*x, 12); 
  b = *y ^ rotr(*y, 12); 
  a = a ^ rotr(a, 3); 
  b = b ^ rotr(b, 3); 
  a = a ^ rotr(*x, 17); 
  b = b ^ rotr(*y, 17); 
  c = a ^ rotr(a, 31); 
  d = b ^ rotr(b, 31); 
  a = a ^ rotr(d, 26); 
  b = b ^ rotr(c, 25); 
  a = a ^ rotr(c, 15); 
  b = b ^ rotr(d, 15); 
  *x = a;
  *y = b;
}

void lbox_layer(uint32_t* state, int r) {
  lbox(&state[0], &state[1+r]);
  lbox(&state[2-r], &state[3]);
}

void add_rc(uint32_t state[LS_ROWS], int round, int shift) {
  for (int i = 0; i < LS_ROWS; i++) {state[i] ^= rc[round][i] << shift;}
}

void dbox_mls_layer(uint32_t state[MLS_BUNDLES][LS_ROWS]) {
  for (int row = 0; row < LS_ROWS; row++) {
    uint32_t x = state[0][row];
    uint32_t y = state[1][row];
    uint32_t z = state[2][row];
    state[0][row] = x ^ y ^ z; 
    state[1][row] = x ^ z; 
    state[2][row] = x ^ y;
  }
}

void shadowenc(uint32_t* x, char* memblock) {

  uint32_t state[MLS_BUNDLES][LS_ROWS]={};

  for (int b = 0; b < MLS_BUNDLES; b++) {
    for (int row = 0; row < LS_ROWS; row++) {
      state[b][row] = x[4*b+row];
    }
  }

  for (int s = 0; s < SHADOW_NS; s++) {
    for (int b = 0; b < MLS_BUNDLES; b++) {
      t_layer(state[b], memblock);
      lbox_layer(state[b], s&1);
      add_rc(state[b], 2 * s, b);
    }
    dbox_mls_layer(state);
    for (int b = 0; b < MLS_BUNDLES; b++) {
      add_rc(state[b], 2 * s + 1, b);
    }
  }

  for (int b = 0; b < MLS_BUNDLES; b++) {
    for (int row = 0; row < LS_ROWS; row++) {
      x[4*b+row] = state[b][row];
    }
  }
}

void generatemessage(uint32_t input[blocknum][SHADOW_NBYTES]){
	for(int i=0; i<blocknum; i++){
    for(int j=0; j<SHADOW_NBYTES;j++){
      input[i][j]=rand();
    }		
	}
}

int main(){

  uint32_t input[blocknum][SHADOW_NBYTES]={};
  uint32_t state[blocknum][SHADOW_NBYTES]={};

  int iter=100000;
  uint64_t start_time, end_time, cyc=0;
  srand(time(0));

  ifstream table("~/table/sbox32.bin", ios::in | ios::binary );

  if(table.is_open()){
    char* memblock=new char[readSize];
    table.read(memblock,readSize);
    table.close();

    generatemessage(input);

    for(int count=0; count<iter; count ++){
      for(int j=0; j<blocknum;j++){
        for(int i=0;i<SHADOW_NBYTES;i++){state[j][i]=input[j][i];}
      }
      
      for(int j=0; j<blocknum;j++){
        start_time=_rdtsc();
        shadowenc(state[j], memblock);
        end_time=_rdtsc();
        cyc+=end_time-start_time;
      }
    }
        
    delete[] memblock;   
  }

  uint64_t aver=cyc/iter;

  cout << "average cost " << aver << " CPU cycles for " << iter << " test" << endl;
	uint64_t cpb = aver/(blocknum*SHADOW_NBYTES*4);
  cout << "average CPB for decryption is " << cpb << endl; 

  return 0;
}