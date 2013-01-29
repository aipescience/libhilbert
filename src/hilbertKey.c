/*  
 *  Copyright (c) 2013, Adrian M. Partl <apartl@aip.de>, 
 *                      eScience team AIP Potsdam
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  See the NOTICE file distributed with this work for additional
 *  information regarding copyright ownership. You may obtain a copy
 *  of the License at
 *
 *  http://www.apache.org/licenses/LICENSE-2.0
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 */

#include "hilbertKey.h"
#include "N10.h"
#include "binaryOps.h"
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

uint64_t getHKeyFromCoord( const int32_t m, const double boxSize, const int32_t dim, const double * point, int * err ) {
	double boxConv = (double)powf(2.0, m) / boxSize;

	//calculate integer values of scaled point to the box coordinate system
	uint64_t iPoint[dim];
	for(int i=0; i < dim; i++) {
		iPoint[i] = (uint64_t)(point[i] * boxConv);
#ifdef VERBOSE
		printf("%i - %llu\n", i, iPoint[i]);
#endif
	}

	return getHKeyFromIntCoord(m, dim, iPoint, err);
}

uint64_t getHKeyFromIntCoord( const int32_t m, const int32_t dim, const uint64_t * point, int * err ) {
	uint64_t result = 0;
	uint64_t TwoPowerOfM = 1 << m;

	uint64_t tmpPoint[dim];
	uint64_t tmp[dim];

	if( dim > HILB_MAX_DIM ) {
		*err = HKEY_ERR_DIM;
		return 0;
	}

	//check data sanity
	for(int i=0; i<dim; i++) {
		assert(point[i] >= 0);
	}

	memcpy(tmpPoint, point, dim * sizeof(uint64_t));
	memset(tmp, 0, dim);

	//clamp larger values to highest possible space on hilbert curve...
	for(int i=0; i<dim; i++) {
		if(tmpPoint[i] > TwoPowerOfM) {
			tmpPoint[i] = TwoPowerOfM - 1;
		}
	}

	//start constructing hilbert key
	for(uint64_t i=0; i<m; i++) {
		//step one: gather upper bits of coordinates
		for(int j=0; j<dim; j++) {
			tmp[j] = IBITS(tmpPoint[j], m-1-i, 1);
#ifdef VERBOSE
		printf("X\n");
		printf("m:%i - %i: %llu, %llu, %llu\n", m, j, i, tmpPoint[j], tmp[j]);
#endif			
		}

		uint64_t upperBitsOfPoint = 0;
		for(int j=0; j<dim; j++) {
			upperBitsOfPoint += tmp[j] << j;
		}

#ifdef VERBOSE	
		printf("X2\n");
		printf("%lli\n", upperBitsOfPoint);
#endif

		//reverse look up value in hilbert template to obtain key number (i.e. H-order)
		uint64_t hOrder = revC[dim-1][upperBitsOfPoint];

#ifdef VERBOSE
		printf("H-Order: %llu\n", hOrder);
#endif		

		//append H-order to result
		result = result << dim;
		result += hOrder;

#ifdef VERBOSE
		printf("R\n");
		printf("%lli\n", result);
#endif

		//perform exchange operation - first figure out which dimensions to exchange (and if...)
		//first isolate the two bits that decode the exchange - do nothing if they are not set
		// x & (-x) (i.e. 01011000 -> 00001000)
		uint32_t * currH = H[dim-1];
		int64_t exchange_gene = currH[hOrder * 2 + 0];
		int64_t reverse_gene = currH[hOrder * 2 + 1];
		uint64_t exDim1 = exchange_gene & (-1 * exchange_gene);
		uint64_t exDim2 = (uint64_t)exchange_gene ^ exDim1;					//this gets the other bit

		//we can have multiple reverses in one go. count them and execute them
		int32_t numReverses = pop64(reverse_gene);
		int64_t tmpReverse_gene = reverse_gene;

		for(int j=0; j<numReverses; j++) {
			int64_t revDim = tmpReverse_gene & ( -1 * tmpReverse_gene );
			tmpReverse_gene = tmpReverse_gene ^ revDim;

#ifdef VERBOSE
			printf("revDim%i: %i %i %i\n", j, tmpReverse_gene, -1 *  tmpReverse_gene, revDim);
#endif

			//reverse the dimensions - first isolating the dimension index by counting 0 bits
			int32_t dimIdx = ntz64(revDim);

			assert(dimIdx < dim);

#ifdef VERBOSE
			printf("Reverse%i %i: %x/%x\n", j, reverse_gene, revDim, dimIdx);
			printf("Reverse%i %i \n", j, tmpPoint[dimIdx]);
#endif

			tmpPoint[dimIdx] = tmpPoint[dimIdx] ^ 0xFFFFFFFFFFFFFFFF;			//reverse operation (i.e. 1011 -> 0100)
		}

		if(exDim1 != 0 && exDim2 != 0) {
			//exchange the dimensions - first isolating the dimension index by counting 0 bits
			int32_t dimIdx1 = ntz64(exDim1);
			int32_t dimIdx2 = ntz64(exDim2);

			assert(dimIdx1 < dim);
			assert(dimIdx2 < dim);

#ifdef VERBOSE
			printf("Exchange: %i, %x/%x <-> %x/%x\n", exchange_gene, exDim1, dimIdx1, exDim2, dimIdx2);
#endif

			uint64_t tmp;

			tmp = tmpPoint[dimIdx1];
			tmpPoint[dimIdx1] = tmpPoint[dimIdx2];
			tmpPoint[dimIdx2] = tmp;
		}

#ifdef VERBOSE	
		printf("C\n");
		printf("The key after iteration %u: ", i);
		for(int j=0; j < dim; j++) {
			printf("dim%i = %llu ", j, tmpPoint[j]);
		}

		printf("result = %llu\n", result);
#endif
	}

#ifdef VERBOSE
	printf("End result = %llu\n", result);
#endif

	*err = HKEY_ERR_OK;
	return result;
}

void getCoordFromHKey( double * outCoord, const int32_t m, const double boxSize, const int32_t dim, const uint64_t key, int * err ) {
	double boxConv = (double)powf(2.0, m) / boxSize;

	//allocate integer array
	uint64_t * result = (uint64_t*)malloc(dim * sizeof(uint64_t));
	if(result == NULL) {
		*err = HKEY_ERR_NOMEM;
		return;
	}
	memset(result, 0, dim * sizeof(uint64_t));

	getIntCoordFromHKey(result, m, dim, key, err);

	if(*err != HKEY_ERR_OK) {
		return;
	}

	for(int i=0; i<dim; i++) {
		outCoord[i] = result[i] / boxConv;
	}

	free(result);

	*err = HKEY_ERR_OK;
	return;
}

void getIntCoordFromHKey( uint64_t * outCoord, const int32_t m, const int32_t dim, const uint64_t key, int * err ) {
	uint64_t tmpKey = key;
	uint64_t flip = 0;

	memset(outCoord, 0, dim * sizeof(uint64_t));

	//checks
	assert(key >= 0);
	assert(key < (pow(2, dim*m)));

	uint64_t lowerNBits = IBITS(tmpKey, 0, dim);
	uint64_t partOnCurve = C[dim-1][lowerNBits];

#ifdef VERBOSE2
		printf("%lli %lli %lli\n", tmpKey, lowerNBits, partOnCurve);
		printf("R\n");
#endif

	for(int i=0; i<dim; i++) {
		outCoord[i] += IBITS(partOnCurve, i, 1);

#ifdef VERBOSE2
		printf("%i %lli\n", i, outCoord[i]);
#endif
	}

	for(int i=1; i<m; i++) {
		lowerNBits = IBITS(tmpKey, dim * i, dim);
		partOnCurve = C[dim-1][lowerNBits];

#ifdef VERBOSE2
		printf("%lli %lli %lli\n", tmpKey, lowerNBits, partOnCurve);
#endif

		//this is needed to perform the bit-flip (adding 111 at the end)
		flip = flip << 1;
		flip = flip + 1;

		//perform exchange operation - first figure out which dimensions to exchange (and if...)
		//first isolate the two bits that decode the exchange - do nothing if they are not set
		// x & (-x) (i.e. 01011000 -> 00001000)
		uint32_t * currH = H[dim-1];
		int64_t exchange_gene = currH[lowerNBits * 2 + 0];
		int64_t reverse_gene = currH[lowerNBits * 2 + 1];
		uint64_t exDim1 = exchange_gene & (-1 * exchange_gene);
		uint64_t exDim2 = (uint64_t)exchange_gene ^ exDim1;					//this gets the other bit

		//we can have multiple reverses in one go. count them and execute them
		int32_t numReverses = pop64(reverse_gene);
		int64_t tmpReverse_gene = reverse_gene;

		if(exDim1 != 0 && exDim2 != 0) {
			//exchange the dimensions - first isolating the dimension index by counting 0 bits
			int32_t dimIdx1 = ntz64(exDim1);
			int32_t dimIdx2 = ntz64(exDim2);

			assert(dimIdx1 < dim);
			assert(dimIdx2 < dim);

			uint64_t tmp;

#ifdef VERBOSE2
			printf("Exchange: %i, %x/%x <-> %x/%x\n", exchange_gene, exDim1, dimIdx1, exDim2, dimIdx2);
#endif

			tmp = outCoord[dimIdx1];
			outCoord[dimIdx1] = outCoord[dimIdx2];
			outCoord[dimIdx2] = tmp;
		}

		for(int j=0; j<numReverses; j++) {
			int64_t revDim = tmpReverse_gene & ( -1 * tmpReverse_gene );
			tmpReverse_gene = tmpReverse_gene ^ revDim;

			//reverse the dimensions - first isolating the dimension index by counting 0 bits
			int32_t dimIdx = ntz64(revDim);

			assert(dimIdx < dim);

#ifdef VERBOSE2
			printf("Reverse%i %i: %x/%x\n", j, reverse_gene, revDim, dimIdx);
			printf("Reverse%i %i \n", j, outCoord[dimIdx]);
#endif

			outCoord[dimIdx] = outCoord[dimIdx] ^ flip;			//reverse operation (i.e. 1011 -> 0100)
		}

#ifdef VERBOSE2
		printf("semi-new outCoord:\n");
		for(int j=0; j<dim; j++) {
			printf("%i %lli\n", j, outCoord[j]);
		}
		printf("R\n");
#endif

		for(int j=0; j<dim; j++) {
			outCoord[j] += IBITS(partOnCurve, j, 1) << i;
#ifdef VERBOSE2
			printf("%i %lli\n", j, outCoord[j]);
#endif			
		}
	}

	*err = HKEY_ERR_OK;
	return;
}