#include <iostream>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <x86intrin.h>

#define blocknum 192 // 3072 bytes

#define readSize  0x400000000
#define CLYDE_NBYTES 4
#define CLYDE_NS 12
#define CLYDE_NR 2 * CLYDE_NS
#define LS_ROWS 4

#define rotr(x,c) (((x) >> (c)) | ((x) << ((32) - (c))))

using namespace std;

static const uint32_t clyde_rc[CLYDE_NR][LS_ROWS] = {
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
};

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

void add_rc(uint32_t state[LS_ROWS], unsigned int round) {
  for (unsigned int i = 0; i < LS_ROWS; i++) {
    state[i] ^= clyde_rc[round][i];
  }
}

void t_layer(uint32_t* state, char* memblock){
  uint64_t tmp;

  for(int i=0; i<LS_ROWS; i++){
    tmp=(uint64_t)4*state[i]; 
    memcpy(&state[i],&memblock[tmp],sizeof(int));
  }
}

void clyde_encrypt(uint32_t* state, char* memblock) {

  for (unsigned int s = 0; s < CLYDE_NS; s++) {
    t_layer(state, memblock);
    lbox_layer(state, s&1);
    add_rc(state, s);
  }
}

void generatemessage(uint32_t input[blocknum][CLYDE_NBYTES]){
	for(int i=0; i<blocknum; i++){
    for(int j=0; j<CLYDE_NBYTES;j++){
      input[i][j]=rand();
    }		
	}
}

int main(){

  uint32_t input[blocknum][CLYDE_NBYTES]={};
  uint32_t state[blocknum][CLYDE_NBYTES]={};

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
      for(int i=0;i<blocknum;i++){
        for(int j=0;j<CLYDE_NBYTES;j++){state[i][j]=input[i][j];}
      }
      
      for(int i=0; i<blocknum;i++){
        start_time=_rdtsc();
        clyde_encrypt(state[i], memblock);
        end_time=_rdtsc();
        cyc+=end_time-start_time;
      }
    }
        
    delete[] memblock;   
  }

  uint64_t aver=cyc/iter;

  cout << "average cost " << aver << " CPU cycles for " << iter << " test" << endl;
	uint64_t cpb = aver/(blocknum*CLYDE_NBYTES*4);
  cout << "average CPB for decryption is " << cpb << endl; 

  return 0;
}
