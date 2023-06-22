#include <iostream>
#include <cstdint>
#include <fstream>
#include <ctime>
#include <x86intrin.h>

#define blocknum 96 // 3072 bytes

#define SATURNIN_CTR_R   8
#define SATURNIN_CTR_D   1

#define round1632 24
#define ROTL8(x,r)  (((x)<<(r)) | (x>>(8-(r))))
#define ROTR8(x,r) (((x)>>(r)) | ((x)<<(8-(r))))
#define newER8(x,y,k_left,k_right) (x^=k_left, y^=k_right, x=ROTR8(x,7), x+=y, y=ROTL8(y,2), y^=x)
#define newDR8(x,y,k_left,k_right) (y^=x, y=ROTR8(y,2), x-=y, x=ROTL8(x,7), y^=k_right, x^=k_left)

using namespace std;

void newSpeck1632Encrypt(uint8_t Pt[], uint8_t Ct[], uint8_t rk[])
{
	uint8_t i;
	Ct[0] = Pt[0]; Ct[1] = Pt[1];
	for (i = 0; i < 2 * round1632; i += 2)
		newER8(Ct[1], Ct[0], rk[i], rk[i + 1]);
}

/*******************   (inv)nonlinear layer of WARX16  ****************************/
void nonlinear(uint16_t arr[16], uint8_t roundkey[]) // need to be modified
{
	uint8_t arrP[16][2],arrC[16][2];
	for (uint16_t i = 0; i < 16; i++)
	{
		arrP[i][0]= (uint8_t)(arr[i] >> 8); //high
		arrP[i][1] = (uint8_t)(arr[i]); //low
		newSpeck1632Encrypt(arrP[i], arrC[i], roundkey);
		arr[i]=(arrC[i][0] << 8) ^ arrC[i][1];
	}  
}

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

void MDS_inv(uint16_t *state)
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

#define MULinv(t0, t1, t2, t3)   do { \
		uint16_t mul_tmp = (t3); \
		(t3) = (t2); \
		(t2) = (t1); \
		(t1) = (t0); \
		(t0) = mul_tmp ^ (t1); \
	} while (0)

	x4 ^= x8; x5 ^= x9; x6 ^= xa; x7 ^= xb; /* B ^= C */
	xc ^= x0; xd ^= x1; xe ^= x2; xf ^= x3; /* D ^= A */
	x8 ^= xc; x9 ^= xd; xa ^= xe; xb ^= xf; /* C ^= D */
	x0 ^= x4; x1 ^= x5; x2 ^= x6; x3 ^= x7; /* A ^= B */
	MULinv(x0, x1, x2, x3);                 /* A = MULinv(A) */
	MULinv(x0, x1, x2, x3);                 /* A = MULinv(A) */
	MULinv(x8, x9, xa, xb);                 /* C = MULinv(C) */
	MULinv(x8, x9, xa, xb);                 /* C = MULinv(C) */
	x4 ^= x8; x5 ^= x9; x6 ^= xa; x7 ^= xb; /* B ^= C */
	xc ^= x0; xd ^= x1; xe ^= x2; xf ^= x3; /* D ^= A */
	MULinv(x4, x5, x6, x7);                 /* B = MULinv(B) */
	MULinv(xc, xd, xe, xf);                 /* D = MULinv(D) */
	x8 ^= xc; x9 ^= xd; xa ^= xe; xb ^= xf; /* C ^= D */
	x0 ^= x4; x1 ^= x5; x2 ^= x6; x3 ^= x7; /* A ^= B */

#undef MULinv

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
		state[ 4 + i] = ((state[ 4 + i] & 0x1111) << 3) | ((state[ 4 + i] & 0xeeee) >> 1);
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

void saturnin_block_decrypt(int R, int D, uint16_t *state)
{
	uint16_t RC0[31]={}, RC1[31]={};
    uint8_t rkeys[48]={0x46,0xb9,0xdd,0x2b,0x0b,0xa8,0x8d,0x13,0x23,0x3b,0x3f,0xeb,0x74,0x3e,0xeb,0x24,0x3f,0xcd,0x52,0xea,0x62,0xb8,0x1b,0x82,0xb5,0x0c,0x27,0x64,0x6e,0xd5,0x76,0x2f,0xd7,0x5d,0xc4,0xdd,0xd8,0xc0,0xf2,0x00,0xcb,0x05,0x01,0x9d,0x67,0xb5,0x92,0xf6};

	int i;

	make_round_constants(R, D, RC0, RC1);

	for (i = R - 1; i >= 0; i --) {

		if ((i & 1) == 0) {
			state[0] ^= RC0[i];
			state[8] ^= RC1[i];
			SR_slice(state);
			MDS_inv(state);
			SR_slice_inv(state);
		} else {
			state[0] ^= RC0[i];
			state[8] ^= RC1[i];
			SR_sheet(state);
			MDS_inv(state);
			SR_sheet_inv(state);
		}
		nonlinear(state,rkeys);
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
        
		for(int j=0;j<blocknum;j++){
			start_time=_rdtsc();
			saturnin_block_decrypt(SATURNIN_CTR_R, SATURNIN_CTR_D, state[j]);
			end_time=_rdtsc();
			cyc+=end_time-start_time;
		}
    }

    uint64_t aver=cyc/iter;

    cout << "average cost " << aver << "CPU cycles for " << iter << " tests \n";
	uint64_t cpb = aver/(2*16*blocknum);
	cout << "average CPB for decryption is " << cpb << "\n";
    return 0;
}