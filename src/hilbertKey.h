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

/*! \file hilbertKey.h
 \brief Hilbert Key generation functions
 
 Hilbert Key algorithm for N-dimensional Hilbert keys. This library implements the
 method described by Chenyang, Hong, Nengchao 2008 IEEE.
 */

#include <stdio.h>
#include <stdint.h>

#ifndef __CLASS_HILBKEY__
#define __CLASS_HILBKEY__

#define HKEY_ERR_DIM   -2 
#define HKEY_ERR_NOMEM -1
#define HKEY_ERR_OK     0

/*! \brief calculate hilbert key from given coordinates in box coordinates (doubles)
 \param const int32_t m:   		hilbert order (max integer dimension: 2**m cells)
 \param const double boxSize:   size of the box for coordinate renormalisation
 \param const int32_t dim:   	number of dimensions
 \param const double * point:   array of size dim with box coordinates of a given point
 \param int * err:   			output variable for error handling
 \return uint64_t hilbert key
 
 Calculates the Hilbert key from coordinates given in box coordinates (doubles). Will
 cast the coordinates to integers according to the hilbert order and calculates hilbert
 key for this "Hilbert cell".*/
uint64_t getHKeyFromCoord( const int32_t m, const double boxSize, const int32_t dim, const double * point, int * err );

/*! \brief calculate hilbert key from given coordinates in coordinates along the hilbert curve
 \param const int32_t m:   		hilbert order (max integer dimension: 2**m cells)
 \param const int32_t dim:   	number of dimensions
 \param const uint64_t * point: array of size dim with coordinates of a given point along hilbert curve (0 < point < 2**m)
 \param int * err:   			output variable for error handling
 \return uint64_t hilbert key
 
 Calculates the Hilbert key from coordinates given in coordinates (int) along the hilbert. 
 curve. If coordinates are larger or smaller 0/2**m, they will clamp to 0/2**m.*/
uint64_t getHKeyFromIntCoord( const int32_t m, const int32_t dim, const uint64_t * point, int * err );

/*! \brief calculate coordinates in box system from a hiven Hilbert key
 \param double * outCoord: 		pre-allocated array for coordinates output
 \param const int32_t m:   		hilbert order (max integer dimension: 2**m cells)
 \param const double boxSize:   size of the box for coordinate renormalisation
 \param const int32_t dim:   	number of dimensions
 \param const uint64_t key: 	hilbert key
 \param int * err:   			output variable for error handling
 \return uint64_t hilbert key
 
 Calculates the coordinates scaled to the box size of a given position along the hilbert curve.
 Result array for the coordinates needs to be allocated before calling this function!*/
void getCoordFromHKey( double * outCoord, const int32_t m, const double boxSize, const int32_t dim, const uint64_t key, int * err );

/*! \brief calculate coordinates from a hiven Hilbert key
 \param double * outCoord: 		pre-allocated array for coordinates output
 \param const int32_t m:   		hilbert order (max integer dimension: 2**m cells)
 \param const int32_t dim:   	number of dimensions
 \param const uint64_t key: 	hilbert key
 \param int * err:   			output variable for error handling
 \return uint64_t hilbert key
 
 Calculates the coordinates of a given position along the hilbert curve.
 Result array for the coordinates needs to be allocated before calling this function!*/
void getIntCoordFromHKey( uint64_t * outCoord, const int32_t m, const int32_t dim, const uint64_t key, int * err );

#endif
