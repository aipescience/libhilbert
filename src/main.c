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

/*! \file main.c
 \brief Test suite for the implementation
 
 Hilbert Key algorithm for N-dimensional Hilbert keys. This library implements the
 method described by Chenyang, Hong, Nengchao 2008 IEEE.
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "hilbertKey.h"
#include <math.h>

void cycleCoordinates_r(const int32_t m, const int32_t dim, const int32_t currDim, uint64_t * point);
void cycleHilbertKey(const int32_t m, const int32_t dim);

int main (int argc, char * const argv[]) {

	if(argc < 4 || argc > 4 + atoi(argv[1])){
	  printf("Usage:\n %s dim m BoxSize [x y z...](optional) \n ", argv[0]);
	  exit(EXIT_SUCCESS);
	}

	int dim = atoi(argv[1]);
	int m = atoi(argv[2]);
	float boxSize = atoi(argv[3]);

	if(argc == 4 + atoi(argv[1])){
		double point[dim];

		for(int i=0; i<dim; i++) {
			point[i] = strtod(argv[4+i], NULL);
			printf("%i: %lf\n", i, point[i]);
		}

		int err;
		uint64_t key = getHKeyFromCoord(m, boxSize, dim, point, &err);

		printf("Hilbert key: %llu\n", key);
	} else {
		uint64_t point[dim];

		for(int i=0; i<pow(2, m); i++) {
			point[0] = (uint64_t) i;
			cycleCoordinates_r(m, dim, 1, point);
		}

		printf("\n");

		cycleHilbertKey(m, dim);
	}

    return EXIT_SUCCESS;
}

void cycleCoordinates_r(const int32_t m, const int32_t dim, const int32_t currDim, uint64_t * point) {
	if(dim == currDim) {
		int err;
		uint64_t key = getHKeyFromIntCoord(m, dim, point, &err);

		printf("%llu ", key);
		for(int i=0; i<dim; i++) {
			printf("%llu ", point[i]);
		}
		printf("\n");
	} else {
		for(int i=0; i<pow(2, m); i++) {
			point[currDim] = (uint64_t) i;
			cycleCoordinates_r(m, dim, currDim + 1, point);
		}
	}
}

void cycleHilbertKey(const int32_t m, const int32_t dim) {
	double tmp = pow(pow(2, m), dim);
	uint64_t numCells = (uint64_t) tmp;
	int err;

	uint64_t * point;
	point = (uint64_t*)malloc(dim * sizeof(uint64_t));

	for(uint64_t i = 0; i<numCells; i++) {
		getIntCoordFromHKey( point, m, dim, i, &err );

		printf("%llu ", i);
		for(int i=0; i<dim; i++) {
			printf("%llu ", point[i]);
		}
		printf("\n");
	}

	free(point);
}