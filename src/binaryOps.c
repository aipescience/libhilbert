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

#include "binaryOps.h"
#include <stdlib.h>
#include <stdio.h>

//number of trailing zeros algorithm using binary search taken from Hacker's delight...
int32_t ntz32(const uint32_t x) {
	uint32_t xCpy = x;
	int32_t n;

	if( xCpy == 0 ) {
		return 32;
	}

	n = 1;
	if(( xCpy & 0x0000FFFF ) == 0) {
		n = n + 16;
		xCpy = xCpy >> 16;
	}
	if(( xCpy & 0x000000FF ) == 0) {
		n = n + 8;
		xCpy = xCpy >> 8;
	}
	if(( xCpy & 0x0000000F ) == 0) {
		n = n + 4;
		xCpy = xCpy >> 4;
	}
	if(( xCpy & 0x00000003 ) == 0) {
		n = n + 2;
		xCpy = xCpy >> 2;
	}

	return n - ( xCpy & 0x00000001 );
}

int32_t ntz64(const uint64_t x) {
	uint64_t xCpy = x;
	int32_t n;

	if( xCpy == 0 ) {
		return 64;
	}

	n = 1;
	if(( xCpy & 0x00000000FFFFFFFF ) == 0) {
		n = n + 32;
		xCpy = xCpy >> 32;
	}
	if(( xCpy & 0x000000000000FFFF ) == 0) {
		n = n + 16;
		xCpy = xCpy >> 16;
	}
	if(( xCpy & 0x00000000000000FF ) == 0) {
		n = n + 8;
		xCpy = xCpy >> 8;
	}
	if(( xCpy & 0x000000000000000F ) == 0) {
		n = n + 4;
		xCpy = xCpy >> 4;
	}
	if(( xCpy & 0x0000000000000003 ) == 0) {
		n = n + 2;
		xCpy = xCpy >> 2;
	}

	return n - ( xCpy & 0x0000000000000001 );
}

int32_t pop32(const uint32_t x) {
	uint32_t xCpy = x;

	xCpy = ( xCpy & 0x55555555 ) + (( xCpy >> 1 )  & 0x55555555 );
	xCpy = ( xCpy & 0x33333333 ) + (( xCpy >> 2 )  & 0x33333333 );
	xCpy = ( xCpy & 0x0F0F0F0F ) + (( xCpy >> 4 )  & 0x0F0F0F0F );
	xCpy = ( xCpy & 0x00FF00FF ) + (( xCpy >> 8 )  & 0x00FF00FF );
	xCpy = ( xCpy & 0x0000FFFF ) + (( xCpy >> 16 ) & 0x0000FFFF );
	return xCpy;
}

int32_t pop64(const uint64_t x) {
	uint64_t xCpy = x;

	xCpy = ( xCpy & 0x5555555555555555 ) + (( xCpy >> 1 )  & 0x5555555555555555 );
	xCpy = ( xCpy & 0x3333333333333333 ) + (( xCpy >> 2 )  & 0x3333333333333333 );
	xCpy = ( xCpy & 0x0F0F0F0F0F0F0F0F ) + (( xCpy >> 4 )  & 0x0F0F0F0F0F0F0F0F );
	xCpy = ( xCpy & 0x00FF00FF00FF00FF ) + (( xCpy >> 8 )  & 0x00FF00FF00FF00FF );
	xCpy = ( xCpy & 0x0000FFFF0000FFFF ) + (( xCpy >> 16 ) & 0x0000FFFF0000FFFF );
	xCpy = ( xCpy & 0x00000000FFFFFFFF ) + (( xCpy >> 32 ) & 0x00000000FFFFFFFF );
	return xCpy;
}
