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

/*! \file binaryOps.h
 \brief Binary number algorithms
 
 Binary number algorithms mostly taken from Hacker's delight by Henry S. Warren, Jr.
 */

#include <stdint.h>

#ifndef __CLASS_BINARYOPS__
#define __CLASS_BINARYOPS__

/*! \brief Implementation of the FORTRAN IBITS function
 \param i:   word
 \param pos: bit position of where to start reading bits
 \param len: length to read in bits
 \return bits
 
 Extracts the len bits at position pos of a word i.*/
#define IBITS(i, pos, len) ((i >> pos) & ~(-1 << len))

/*! \brief Number of trailing zeros 32bits
 \param const uint32_t x:   wird
 \return int32_t number of trailing zeros
 
 Counts the number of trailing zeros of a 32 bit unsigned integer.*/
int32_t ntz32(const uint32_t x);

/*! \brief Number of trailing zeros 64bits
 \param const uint64_t x:   word
 \return int32_t number of trailing zeros
 
 Counts the number of trailing zeros of a 64 bit unsigned integer.*/
int32_t ntz64(const uint64_t x);


/*! \brief Number of 1-bits in a given 32 bit word
 \param const uint32_t x:   variable
 \return int32_t number of 1-bits
 
 Counts the number of 1-bits in a 32 bit unsigned integer.*/
int32_t pop32(const uint32_t x);

/*! \brief Number of 1-bits in a given 64 bit word
 \param const uint64_t x:   variable
 \return int32_t number of 1-bits
 
 Counts the number of 1-bits in a 64 bit unsigned integer.*/
int32_t pop64(const uint64_t x);

#endif
