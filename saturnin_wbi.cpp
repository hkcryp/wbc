#include <iostream>
#include <cstdint>
#include <fstream>
#include <ctime>
#include <x86intrin.h>
#include "tbox.h"

#define blocknum 96

#define SATURNIN_CTR_R   8
#define SATURNIN_CTR_D   1

using namespace std;

void make_round_constants(int R, int D, uint16_t *RC0, uint16_t *RC1)
{
	uint16_t x0, x1;
	int n;

	x0 = x1 = D + (R << 4) + 0xFE00;

	for (n = 0; n < R; n ++) {
		int i;

		for (i = 0; i < 16; i ++) {
			x0 = (x0 << 1) ^ (0x2D & -(x0 >> 15));
			x1 = (x1 << 1) ^ (0x53 & -(x1 >> 15));
		}
		RC0[n] = x0;
		RC1[n] = x1;
	}
}

void MDS(uint16_t *state)
{
	uint16_t x0, x1, x2, x3, x4, x5, x6, x7;
	uint16_t x8, x9, xa, xb, xc, xd, xe, xf;

	x0 = state[0x0];
	x1 = state[0x1];
	x2 = state[0x2];
	x3 = state[0x3];
	x4 = state[0x4];
	x5 = state[0x5];
	x6 = state[0x6];
	x7 = state[0x7];
	x8 = state[0x8];
	x9 = state[0x9];
	xa = state[0xa];
	xb = state[0xb];
	xc = state[0xc];
	xd = state[0xd];
	xe = state[0xe];
	xf = state[0xf];

#define MUL(t0, t1, t2, t3)   do { \
		uint16_t mul_tmp = (t0); \
		(t0) = (t1); \
		(t1) = (t2); \
		(t2) = (t3); \
		(t3) = mul_tmp ^ (t0); \
	} while (0)

	x8 ^= xc; x9 ^= xd; xa ^= xe; xb ^= xf;  /* C ^= D */
	x0 ^= x4; x1 ^= x5; x2 ^= x6; x3 ^= x7;  /* A ^= B */
	MUL(x4, x5, x6, x7);                     /* B = MUL(B) */
	MUL(xc, xd, xe, xf);                     /* D = MUL(D) */
	x4 ^= x8; x5 ^= x9; x6 ^= xa; x7 ^= xb;  /* B ^= C */
	xc ^= x0; xd ^= x1; xe ^= x2; xf ^= x3;  /* D ^= A */
	MUL(x0, x1, x2, x3);                     /* A = MUL(A) */
	MUL(x0, x1, x2, x3);                     /* A = MUL(A) */
	MUL(x8, x9, xa, xb);                     /* C = MUL(C) */
	MUL(x8, x9, xa, xb);                     /* C = MUL(C) */
	x8 ^= xc; x9 ^= xd; xa ^= xe; xb ^= xf;  /* C ^= D */
	x0 ^= x4; x1 ^= x5; x2 ^= x6; x3 ^= x7;  /* A ^= B */
	x4 ^= x8; x5 ^= x9; x6 ^= xa; x7 ^= xb;  /* B ^= C */
	xc ^= x0; xd ^= x1; xe ^= x2; xf ^= x3;  /* D ^= A */

#undef MUL

	state[0x0] = x0;
	state[0x1] = x1;
	state[0x2] = x2;
	state[0x3] = x3;
	state[0x4] = x4;
	state[0x5] = x5;
	state[0x6] = x6;
	state[0x7] = x7;
	state[0x8] = x8;
	state[0x9] = x9;
	state[0xa] = xa;
	state[0xb] = xb;
	state[0xc] = xc;
	state[0xd] = xd;
	state[0xe] = xe;
	state[0xf] = xf;
}

void SR_slice(uint16_t *state)
{
	int i;

	for (i = 0; i < 4; i ++) {
		state[ 4 + i] = ((state[ 4 + i] & 0x7777) << 1) | ((state[ 4 + i] & 0x8888) >> 3);
		state[ 8 + i] = ((state[ 8 + i] & 0x3333) << 2) | ((state[ 8 + i] & 0xcccc) >> 2);
		state[12 + i] = ((state[12 + i] & 0x1111) << 3) | ((state[12 + i] & 0xeeee) >> 1);
	}
}

void SR_slice_inv(uint16_t *state)
{
	int i;

	for (i = 0; i < 4; i ++) {
		state[ 4 + i] = ((state[ 4 + i] & 0x1111) << 3)	| ((state[ 4 + i] & 0xeeee) >> 1);
		state[ 8 + i] = ((state[ 8 + i] & 0x3333) << 2)	| ((state[ 8 + i] & 0xcccc) >> 2);
		state[12 + i] = ((state[12 + i] & 0x7777) << 1)	| ((state[12 + i] & 0x8888) >> 3); 
	}
}

void SR_sheet(uint16_t *state)
{
	int i;

	for (i = 0; i < 4; i ++) {
		state[ 4 + i] = ((state[ 4 + i] <<  4) | (state[ 4 + i] >> 12));
		state[ 8 + i] = ((state[ 8 + i] <<  8) | (state[ 8 + i] >>  8));
		state[12 + i] = ((state[12 + i] << 12) | (state[12 + i] >>  4));
	}
}

void SR_sheet_inv(uint16_t *state)
{
	int i;

	for (i = 0; i < 4; i ++) {
		state[ 4 + i] = ((state[ 4 + i] << 12) | (state[ 4 + i] >>  4));
		state[ 8 + i] = ((state[ 8 + i] <<  8) | (state[ 8 + i] >>  8)); 
		state[12 + i] = ((state[12 + i] <<  4) | (state[12 + i] >> 12)); 
	}
}

void t_layer(uint16_t *state){
    for(int i=0;i<16;i++){state[i]=invT[state[i]];}
}

void saturnin_block_encrypt(int R, int D, uint16_t *xb)
{
	uint16_t RC0[31], RC1[31];

	int i;

	make_round_constants(R, D, RC0, RC1);

	for (i = 0; i < R; i ++) {

		t_layer(xb);

		if ((i & 1) == 0) {
			SR_slice(xb);
			MDS(xb);
			SR_slice_inv(xb);
			xb[0] ^= RC0[i]; 
			xb[8] ^= RC1[i]; 
		} else {
			SR_sheet(xb);
			MDS(xb);
			SR_sheet_inv(xb);
			xb[0] ^= RC0[i]; 
			xb[8] ^= RC1[i]; 
		}
	}
}

void generatemessage(uint16_t input[blocknum][16]){
	for(int i=0; i<blocknum; i++){
		for(int j=0; j<16; j++){input[i][j]=rand();}
	}
}

int main(){
    uint16_t input[blocknum][16] = {};
	uint16_t state[blocknum][16] = {};

    int iter=100000;
	srand(time(0));

    uint64_t start_time, end_time, cyc=0;

	generatemessage(input);

    for(int count=0; count<iter;count++){
        for(int i=0;i<blocknum;i++){
			for(int j=0;j<16;j++){state[i][j]=input[i][j];}
		}
        start_time=_rdtsc();
		for(int j=0;j<blocknum;j++){
			saturnin_block_encrypt(SATURNIN_CTR_R, SATURNIN_CTR_D, state[j]);
		}
        end_time=_rdtsc();
        cyc+=end_time-start_time;
    }

    uint64_t aver=cyc/iter;

    cout << "average cost " << aver << "CPU cycles for " << iter << " tests \n";
	uint64_t cpb = aver/(2*16*blocknum);
	cout << "average CPB for decryption is " << cpb << "\n";
    return 0;
}