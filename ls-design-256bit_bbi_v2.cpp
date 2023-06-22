#include <iostream>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <fstream>
#include <x86intrin.h>

#define blocknum 96

#define SHADOW_NBYTES 8
#define SHADOW_NS 14 
#define TABLE_NS 14
#define LS_ROWS 4   
#define LS_ROW_BYTES 4 
#define MLS_BUNDLES (SHADOW_NBYTES / LS_ROWS)

#define xtime(a) (((a) << (1)) ^ (((a) >> (3)) * (0x13)))
#define rotr(x,c) (((x) >> (c)) | ((x) << ((32) - (c))))

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

uint32_t shadow_rc[SHADOW_NS][LS_ROWS] = {
  { 1, 0, 0, 1 },                                       // 0 { 1, 0, 0, 1 }
  { 0, 2, 0, 0 },                                       // 1 { 0, 1, 0, 0 }
  { 0, 0, 4, 0 },                                       // 2 { 0, 0, 1, 0 }
  { 0, 0, 0, 8 },                                       // 3 { 0, 0, 0, 1 }
  { 0x10, 0, 0, 0 },                                    // 4 { 1, 0, 0, 0 }
  { 0x20, 0x20, 0, 0 },                                 // 5 { 1, 1, 0, 0 }
  { 0x40, 0x40, 0x40, 0 },                              // 6 { 1, 1, 1, 0 }
  { 0x80, 0x80, 0x80, 0x80 },                           // 7 { 1, 1, 1, 1 }
  { 0, 0x100, 0x100, 0x100 },                           // 8 { 0, 1, 1, 1 }
  { 0x200, 0, 0x200, 0x200 },                           // 9 { 1, 0, 1, 1 }
  { 0, 0x400, 0, 0x400 },                               // 10 { 0, 1, 0, 1 }
  { 0x800, 0, 0x800, 0 },                               // 11 { 1, 0, 1, 0 }
  { 0x1000, 0x1000, 0, 0x1000 },                        // 12 { 1, 1, 0, 1 }
  { 0, 0x2000, 0x2000, 0 },                             // 13 { 0, 1, 1, 0 }
  // { 0, 0, 0x4000, 0x4000 },                             // 14 { 0, 0, 1, 1 }
};

uint32_t shadow_drc[SHADOW_NS][2*LS_ROWS] = {
  {0, 0, 1, 0, 1, 0, 1, 1},                              // 0 {0, 0, 1, 0, 1, 0, 1, 1}
  {0, 0, 0, 2, 0, 2, 0, 2},                              // 1 {0, 0, 0, 1, 0, 1, 0, 1}
  {4, 0, 0, 0, 4, 0, 4, 0},                              // 2 {1, 0, 0, 0, 1, 0, 1, 0}
  {0, 8, 0, 0, 0, 8, 0, 8},                              // 3 {0, 1, 0, 0, 0, 1, 0, 1}
  {0x10, 0, 0x10, 0, 0, 0, 0x10, 0},                     // 4 {1, 0, 1, 0, 0, 0, 1, 0}
  {0, 0x20, 0, 0x20, 0, 0, 0, 0x20},                     // 5 {0, 1, 0, 1, 0, 0, 0, 1}
  {0x40, 0, 0x40, 0, 0x40, 0, 0, 0},                     // 6 {1, 0, 1, 0, 1, 0, 0, 0}
  {0, 0x80, 0, 0x80, 0, 0x80, 0, 0},                     // 7 {0, 1, 0, 1, 0, 1, 0, 0}
  {0, 0, 0x100, 0, 0x100, 0, 0x100, 0},                  // 8 {0, 0, 1, 0, 1, 0, 1, 0}
  {0x200, 0, 0, 0x200, 0, 0x200, 0, 0x200},              // 9 {1, 0, 0, 1, 0, 1, 0, 1}
  {0x400, 0x400, 0, 0, 0x400, 0, 0x400, 0},              // 10 {1, 1, 0, 0, 1, 0, 1, 0}
  {0, 0x800, 0x800, 0, 0, 0x800, 0, 0x800},              // 11 {0, 1, 1, 0, 0, 1, 0, 1}
  {0, 0, 0x1000, 0x1000, 0, 0, 0x1000, 0},               // 12 {0, 0, 1, 1, 0, 0, 1, 0}
  {0, 0, 0, 0x2000, 0x2000, 0, 0, 0x2000},               // 13 {0, 0, 0, 1, 1, 0, 0, 1}
  // {0, 0, 0, 0, 0x4000, 0x4000, 0, 0},                    // 14 {0, 0, 0, 0, 1, 1, 0, 0}
};

uint8_t shadow_trc[TABLE_NS][8] = {
{0, 1, 1, 1, 0, 0, 0, 1}, //0
{0, 0, 1, 0, 1, 1, 1, 0}, //1
{1, 1, 0, 0, 0, 0, 1, 1}, //2
{0, 0, 0, 0, 0, 1, 1, 0}, //3
{1, 1, 0, 1, 0, 1, 0, 1}, //4
{0, 0, 0, 1, 0, 0, 1, 0}, //5
{1, 1, 1, 1, 0, 0, 1, 0}, //6
{1, 1, 0, 0, 0, 1, 0, 0}, //7
{0, 0, 0, 0, 1, 0, 1, 0}, //8
{1, 1, 1, 1, 1, 1, 1, 1}, //9
{0, 0, 1, 1, 0, 1, 1, 1}, //10
{0, 0, 0, 1, 0, 1, 1, 1}, //11
{1, 1, 0, 0, 1, 1, 0, 1}, //12
{0, 0, 0, 1, 1, 1, 1, 1}, //13
// {0, 0, 0, 0, 0, 0, 0, 1}, //14
// {1, 1, 0, 1, 1, 0, 0, 1}, //15
// {0, 0, 1, 1, 1, 0, 0, 0}, //16
// {1, 1, 0, 1, 0, 1, 1, 1}, //17
// {0, 0, 1, 0, 0, 0, 0, 1}, //18
// {0, 0, 0, 0, 0, 0, 1, 1}, //19
// {1, 1, 1, 0, 1, 0, 1, 0}, //20
// {1, 1, 0, 0, 1, 0, 0, 1}, //21
// {1, 1, 1, 1, 1, 0, 0, 1}, //22
// {1, 1, 1, 0, 0, 0, 1, 0}, //23
// {0, 0, 0, 0, 0, 1, 0, 1}, //24
// {0, 0, 1, 1, 1, 1, 1, 1}, //25
// {1, 1, 0, 1, 1, 0, 1, 1}, //26
// {0, 0, 0, 0, 1, 0, 1, 1}, //27
// {0, 0, 1, 0, 0, 1, 1, 0}, //28
// {0, 0, 0, 0, 1, 1, 1, 1}, //29
// {1, 1, 0, 0, 0, 0, 0, 0}, //30
// {1, 1, 1, 0, 1, 1, 0, 0}, //31
};

uint8_t invSBox[256] = {0x8C, 0x21, 0xDB, 0xB3, 0x34, 0x32, 0xFE, 0xA7, 0x53, 0x7B, 0xC7, 0x6B, 0xA9, 0xA5, 0x22, 0x50,
                        0x97, 0xA3, 0xAA, 0x31, 0x16, 0xC9, 0x9F, 0x15, 0xF2, 0xE9, 0xDC, 0xF9, 0xE3, 0x84, 0xD9, 0x4A,
                        0x00, 0x6D, 0x3C, 0x0F, 0x38, 0xB7, 0x82, 0x4B, 0xAB, 0x6E, 0x6F, 0x87, 0x41, 0x69, 0x5E, 0x99,
                        0x12, 0x36, 0x86, 0x04, 0x83, 0x05, 0x88, 0x30, 0x39, 0x25, 0x55, 0x9C, 0x7A, 0x23, 0xBC, 0x62,
                        0x2D, 0x74, 0xE8, 0x46, 0x61, 0xC3, 0xAD, 0x42, 0xA0, 0x56, 0x26, 0x1E, 0xAC, 0xE6, 0xD3, 0x5D,
                        0x76, 0x0E, 0x09, 0xFC, 0x3B, 0x78, 0x5C, 0x48, 0x91, 0xEC, 0x7D, 0x64, 0x4E, 0x57, 0x68, 0x2F,
                        0x45, 0xCC, 0x63, 0x3E, 0x89, 0x5A, 0xC5, 0x6A, 0x2C, 0x5F, 0x0A, 0x66, 0x20, 0x9E, 0x2B, 0x28,
                        0x7F, 0xD7, 0xE1, 0x75, 0x72, 0x40, 0x77, 0x51, 0x96, 0x54, 0x08, 0x3D, 0x5B, 0x7C, 0x71, 0x93,
                        0xB8, 0x85, 0x35, 0x27, 0x80, 0x1C, 0x2A, 0x33, 0x65, 0x37, 0xC1, 0xAF, 0xEF, 0x01, 0xA6, 0x8E,
                        0x59, 0xBF, 0x7E, 0xDD, 0x9A, 0xED, 0x11, 0x79, 0x2E, 0xCD, 0xF0, 0x95, 0x9D, 0x3A, 0x17, 0x6C,
                        0xB4, 0x49, 0x10, 0xCB, 0x0C, 0xFB, 0x06, 0x8F, 0x0D, 0xB0, 0x29, 0x13, 0x47, 0x4D, 0x8A, 0xF5,
                        0xCE, 0xA8, 0x02, 0xBA, 0xFD, 0xA1, 0x24, 0xEE, 0xF7, 0x81, 0xBB, 0xB2, 0xF6, 0x3F, 0x90, 0xD4,
                        0x8B, 0xC2, 0x44, 0xC0, 0x67, 0xDF, 0x0B, 0xF4, 0x14, 0xC8, 0xA2, 0xD8, 0x98, 0x60, 0xE5, 0xB1,
                        0xFA, 0xD0, 0x4F, 0xD2, 0xD5, 0xBE, 0x70, 0xD6, 0x1F, 0xCA, 0x03, 0xDA, 0x92, 0x1B, 0xC4, 0xEB,
                        0x73, 0xEA, 0x1D, 0xF8, 0xCF, 0xE4, 0xF3, 0x4C, 0x18, 0x43, 0xDE, 0xE0, 0x94, 0x58, 0x8D, 0xB6,
                        0xF1, 0x9B, 0xE7, 0x19, 0xAE, 0xC6, 0xB9, 0xBD, 0x1A, 0xE2, 0xA4, 0xD1, 0xB5, 0x52, 0xFF, 0x07};

uint16_t xTimeTable[3][16] = {};
uint32_t lBoxTable[8][256] = {};

void GenerateXTimeTable()
{
	for (int i = 0; i < 16; i++)
	{
		xTimeTable[0][i] = xtime(i);
		xTimeTable[1][i] = xtime(xtime(i));
		xTimeTable[2][i] = xtime(xtime(xtime(i)));
	}
}

void GenerateLBoxTable()
{
	for (int i = 0; i < 256; i++)
	{
		for (int j = 0; j < 8; j++)
		{
			lBoxTable[0][i] |= (i & (1 << j)) << (j * 3);
		}
	}

	for (int i = 0; i < 256; i++)
	{
		for (int j = 0; j < 8; j++)
		{
			lBoxTable[1][i] |= (i & (1 << j)) << ((j * 3) + 1);
		}
	}

	for (int i = 0; i < 256; i++)
	{
		for (int j = 0; j < 8; j++)
		{
			lBoxTable[2][i] |= (i & (1 << j)) << ((j * 3) + 2);
		}
	}

	for (int i = 0; i < 256; i++)
	{
		for (int j = 0; j < 8; j++)
		{
			lBoxTable[3][i] |= (i & (1 << j)) << ((j * 3) + 3);
		}
	}

	for (int i = 0; i < 256; i++)
	{	
		lBoxTable[4][i] |= (i & 1) | ((i & 0x2) << 7) | ((i & 0x4) << 14) | ((i & 0x8) << 21);
		lBoxTable[4][i] |= ((i & 0x10) >> 3) | ((i & 0x20) << 4) | ((i & 0x40) << 11) | ((i & 0x80) << 18);
	}

	for (int i = 0; i < 256; i++)
	{
		lBoxTable[5][i] |= ((i & 1) << 2) | ((i & 0x2) << 9) | ((i & 0x4) << 16) | ((i & 0x8) << 23);
		lBoxTable[5][i] |= ((i & 0x10) >> 1) | ((i & 0x20) << 6) | ((i & 0x40) << 13) | ((i & 0x80) << 20);
	}

	for (int i = 0; i < 256; i++)
	{
		lBoxTable[6][i] |= ((i & 1) << 4) | ((i & 0x2) << 11) | ((i & 0x4) << 18) | ((i & 0x8) << 25);
		lBoxTable[6][i] |= ((i & 0x10) << 1) | ((i & 0x20) << 8 ) | ((i & 0x40) << 15) | ((i & 0x80) << 22);
	}

	for (int i = 0; i < 256; i++)
	{
		lBoxTable[7][i] |= ((i & 1) << 6) | ((i & 0x2) << 13) | ((i & 0x4) << 20) | ((i & 0x8) << 27);
		lBoxTable[7][i] |= ((i & 0x10) << 3) | ((i & 0x20) << 10) | ((i & 0x40) << 17) | ((i & 0x80) << 24);
	}
}
 
void table_lbox_inv(uint8_t *state){

	uint8_t l[8] = {};
  uint8_t s[8] = {};
	
  for(int i=0; i<4; i++){
    l[2*i] = state[i] & 0xF;
    l[2*i + 1] = (state[i] & 0xF0) >> 4;
  }

    s[7] = l[7] ^ l[5] ^ l[4] ^ l[3] ^ l[2] ^ l[1] ^ xTimeTable[0][l[5] ^ l[4] ^ l[2] ^ l[1]] ^ xTimeTable[1][l[4] ^ l[3] ^ l[2]] ^ xTimeTable[2][l[6] ^ l[4] ^ l[2] ^ l[0]];
    s[6] = l[6] ^ l[5] ^ l[3] ^ l[0] ^ xTimeTable[0][l[5] ^ l[4] ^ l[0]] ^ xTimeTable[1][l[6] ^ l[2] ^ l[1] ^ l[0]] ^ xTimeTable[2][l[7] ^ l[6] ^ l[0]];
    s[5] = l[7] ^ l[6] ^ l[5] ^ l[4] ^ l[2] ^ l[0] ^ xTimeTable[0][l[7] ^ l[5] ^ l[2] ^ l[1]] ^ xTimeTable[1][l[7] ^ l[5] ^ l[3] ^ l[1] ^ l[0]] ^ xTimeTable[2][l[7] ^ l[6] ^ l[5] ^ l[4] ^ l[2]];
    s[4] = l[7] ^ l[6] ^ l[4] ^ l[3] ^ xTimeTable[0][l[5] ^ l[3] ^ l[2]] ^ xTimeTable[1][l[7] ^ l[5] ^ l[1]] ^ xTimeTable[2][l[4] ^ l[3] ^ l[0]];
    s[3] = l[6] ^ l[4] ^ l[3] ^ l[1] ^ xTimeTable[0][l[5] ^ l[4] ^ l[3] ^ l[2]] ^ xTimeTable[1][l[4] ^ l[3]] ^ xTimeTable[2][l[7] ^ l[6] ^ l[5] ^ l[2] ^ l[1] ^ l[0]];
    s[2] = l[4] ^ l[3] ^ l[1] ^ l[0] ^ xTimeTable[0][l[5] ^ l[4] ^ l[2]] ^ xTimeTable[1][l[6] ^ l[2] ^ l[0]] ^ xTimeTable[2][l[7] ^ l[4] ^ l[3]];
    s[1] = l[7] ^ l[5] ^ l[3] ^ l[2] ^ l[1] ^ l[0] ^ xTimeTable[0][l[6] ^ l[5] ^ l[2] ^ l[0]] ^ xTimeTable[1][l[7] ^ l[6] ^ l[4] ^ l[2] ^ l[0]] ^ xTimeTable[2][l[5] ^ l[3] ^ l[2] ^ l[1] ^ l[0]];
    s[0] = l[7] ^ l[4] ^ l[2] ^ l[1] ^ xTimeTable[0][l[7] ^ l[3] ^ l[2]] ^ xTimeTable[1][l[7] ^ l[6] ^ l[5] ^ l[1]] ^ xTimeTable[2][l[7] ^ l[1] ^ l[0]];


  for(int i=0; i<4; i++){
    state[i] = (s[2*i+1] << 4) | s[2*i];
  }
}

void table_sbox_layer_inv(uint8_t* state) {
  for(int i=0; i<4; i++){state[i] = invSBox[state[i]];}
}

void addroundkey_layer(uint8_t* state, uint8_t* roundkey) { 
  for(int i=0; i<4; i++){state[i] ^= roundkey[i];}
}

void table_add_rc(uint8_t* state, int round) {

  for(int i=0; i<4; i++){state[i] ^= shadow_trc[round][i];}
	
}

void invt_layer(uint32_t* input){

  uint8_t key[0x20][4] = {{0x46, 0xb9, 0xdd, 0x2b}, {0x0b, 0xa8, 0x8d, 0x13}, {0x23, 0x3b, 0x3f, 0xeb}, {0x74, 0x3e, 0xeb, 0x24}, {0x3f, 0xcd, 0x52, 0xea}, {0x62, 0xb8, 0x1b, 0x82}, {0xb5, 0x0c, 0x27, 0x64}, {0x6e, 0xd5, 0x76, 0x2f}, 
                          {0xd7, 0x5d, 0xc4, 0xdd}, {0xd8, 0xc0, 0xf2, 0x00}, {0xcb, 0x05, 0x01, 0x9d}, {0x67, 0xb5, 0x92, 0xf6}, {0xfc, 0x82, 0x1c, 0x49}, {0x47, 0x9a, 0xb4, 0x86}, {0x40, 0x29, 0x2e, 0xac}, {0xb3, 0xb7, 0xc4, 0xbe}, 
                          {0x46, 0xb9, 0xdd, 0x2b}, {0x0b, 0xa8, 0x8d, 0x13}, {0x23, 0x3b, 0x3f, 0xeb}, {0x74, 0x3e, 0xeb, 0x24}, {0x3f, 0xcd, 0x52, 0xea}, {0x62, 0xb8, 0x1b, 0x82}, {0xb5, 0x0c, 0x27, 0x64}, {0x6e, 0xd5, 0x76, 0x2f},
                          {0xd7, 0x5d, 0xc4, 0xdd}, {0xd8, 0xc0, 0xf2, 0x00}, {0xcb, 0x05, 0x01, 0x9d}, {0x67, 0xb5, 0x92, 0xf6}, {0xfc, 0x82, 0x1c, 0x49}, {0x47, 0x9a, 0xb4, 0x86}, {0x40, 0x29, 0x2e, 0xac}, {0xb3, 0xb7, 0xc4, 0xbe}};

      
    uint8_t state[4] = {};
    uint32_t st = 0;

    for(int j=0; j<LS_ROWS; j++){

      st = lBoxTable[0][input[j] & 0xFF] | lBoxTable[1][(input[j] & 0x0000FF00) >> 8] | lBoxTable[2][(input[j] & 0x00FF0000) >> 16] | lBoxTable[3][(input[j] & 0xFF000000) >> 24];

      state[0] = (st & 0xFF) ;
      state[1] = ((st & 0x0000FF00) >> 8) ;
      state[2] = ((st & 0x00FF0000) >> 16);
      state[3] = ((st & 0xFF000000) >> 24);

      for (int s = (TABLE_NS -1); s >= 0; s--) {
        table_add_rc(state, s);
        table_lbox_inv(state);
        table_sbox_layer_inv(state);
        addroundkey_layer(state, key[s]);     
      }

      st = (state[3] << 24) | (state[2] << 16) | (state[1] << 8) | (state[0]);
      input[j] = lBoxTable[4][st & 0xFF] | lBoxTable[5][(st & 0x0000FF00) >> 8] | lBoxTable[6][(st & 0x00FF0000) >> 16] | lBoxTable[7][(st & 0xFF000000) >> 24];
    }
}

static void lbox_inv(uint32_t* x, uint32_t* y) {
  uint32_t a, b, c, d;
  a = *x ^ rotr(*x, 25);
  b = *y ^ rotr(*y, 25);
  c = *x ^ rotr(a, 31);
  d = *y ^ rotr(b, 31);
  c = c ^ rotr(a, 20);
  d = d ^ rotr(b, 20);
  a = c ^ rotr(c, 31);
  b = d ^ rotr(d, 31);
  c = c ^ rotr(b, 26);
  d = d ^ rotr(a, 25);
  a = a ^ rotr(c, 17);
  b = b ^ rotr(d, 17);
  a = rotr(a, 16);
  b = rotr(b, 16);
  *x = a;
  *y = b;
}

static void lbox_layer_inv(uint32_t* state, int r) {
    lbox_inv(&state[0], &state[1+r]);
    lbox_inv(&state[2-r], &state[3]);
}

void add_rc(uint32_t state[LS_ROWS], int round) {
  for (int i = 0; i < LS_ROWS; i++) { state[i] ^= shadow_rc[round][i]; }
}

void add_drc(uint32_t state[LS_ROWS], int round) {
  for (int i = 0; i < LS_ROWS; i++) { state[i] ^= shadow_drc[round][i]; }
}

void dbox_mls_layer_inv(uint32_t* state, uint32_t* s) {

  for(int i=0; i<4; i++){
    s[i] ^= rotr(state[i],19);
    state[i] ^= rotr(s[i],15);
  }
}

void shadowdec(uint32_t* x) {

  uint32_t state[MLS_BUNDLES][LS_ROWS]={};

  for (int b = 0; b < MLS_BUNDLES; b++) {
    for (int row = 0; row < LS_ROWS; row++) {
      state[b][row] = x[4*b+row];
    }
  }

  for (int s = SHADOW_NS-1; s >= 0 ; s--) {
    for (int b = MLS_BUNDLES-1; b >= 0; b--) {
      add_drc(state[b], s);
    }
    dbox_mls_layer_inv(state[0],state[1]);
    for (int b = MLS_BUNDLES-1; b >= 0; b--) {
      add_rc(state[b], s);
      lbox_layer_inv(state[b], s&1);
      invt_layer(state[b]);     
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

  generatemessage(input);
  GenerateXTimeTable();
	GenerateLBoxTable();

  for(int count=0; count<iter; count ++){
    for(int j=0; j<blocknum;j++){
      for(int i=0;i<SHADOW_NBYTES;i++){state[j][i]=input[j][i];}
    }    
    start_time=_rdtsc();
    for(int j=0;j<blocknum;j++){
      shadowdec(state[j]);
      }    
    end_time=_rdtsc();
    cyc+=end_time-start_time;
  }

  uint64_t aver=cyc/iter;

  cout << "average cost " << aver << " CPU cycles for " << iter << " test" << endl;
	uint64_t cpb = aver/(blocknum*SHADOW_NBYTES*4);
  cout << "average CPB for decryption is " << cpb << endl; 

  return 0;
}

