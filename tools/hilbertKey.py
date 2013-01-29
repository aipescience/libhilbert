""" 
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
"""

"""
 * Script for generating hilbert generation genes. This is used for semi-auto-generating
 * the C header files. Will generate numbers up to given dimension.
 * WARNING: N larger than 20 will produce HUGE (multiple MB) sized header files...
 *
 * Usage: python hilbertKey dim (number of dimensions)
"""

import sys
import numpy as np
import math as m

def calC_r(N):
	result = np.zeros(shape=(2**N, N), dtype=np.uint8)

	if(N == 1):
		result[0][0] = 0
		result[1][0] = 1
		return result

	C_n1 = calC_r(N-1)

	for i in range(0, 2**N / 2):
		result[i][0] = 0
		result[i][1:] = C_n1[i]

	for i in range(0, 2**N / 2):
		result[2**N / 2 + i][0] = 1
		result[2**N / 2 + i][1:] = C_n1[2**N / 2 - i - 1]

	return result


def calH1_N(i, C, N, H):
	result = np.zeros(shape=(2, N), dtype=np.uint8)

	if ( i < 2**N / 2 ):
		if ( i == 0 ):
			result[0] = C[0]

			result[1] = result[0] ^ (C[0] ^ C[1])
		else:
			result[0] = H[i-1][1] ^ (C[i-1] ^ C[i])
			tmp = result[0] & (C[i] ^ C[i+1])
			tmp2 = C[i] & (C[i] ^ C[i+1])
			if(np.all(tmp == tmp2)):
				result[1] = result[0] ^ (C[i] ^ C[i+1])
			else:
				#we have multiple options here, just choose one
				if 1==0:
					print "Special...", i
					print tmp2
					print tmp
					print H[i-1][1]
					print "1 ", (C[i] ^ C[i+1])
					print C[i]
					print C[i+1]
					print "R0 ",result[0]

				#build array of all options and prune
				options = np.zeros(shape=(N, N), dtype=np.uint8)

				count = 0
				for j in range(N):
					tmpArray = np.zeros(shape=(N), dtype=np.uint8)
					tmpArray[j] = 1
					if(np.all(tmpArray == (C[i] ^ C[i+1]))):
						continue
					else:
						options[count] = result[0] ^ tmpArray
						count += 1

				options = options[:count]
				#print options
				result[1] = options[0]

	else:
		if( i == 2**N / 2 ):
			tmp = np.zeros(shape=(N), dtype=np.uint8)
			tmp[0] = 1  #original 0x000000001<<n-1
			result[0] = H[2**N-i-1][1] ^ tmp
			result[1] = result[0] ^ (H[2**N-i-1][0] ^ H[2**N-i-1][1])
		else:
			result[0] = H[i-1][1] ^ (H[2**N-i][0] ^ H[2**N-i-1][1])
			result[1] = result[0] ^ (H[2**N-i-1][0] ^ H[2**N-i-1][1])

	#print "Res; ", result

	return result

def calcH(C, N):
	H = np.zeros(shape=(2**N, 2, N), dtype=np.uint8)
	for i in range(2**N):
		H[i] = calH1_N(i, C, N, H)

	return H

def calcG(C, H, N):
	G = np.zeros(shape=(2**N, 2, N), dtype=np.uint8)
	for i in range(2**N):
		G[i][0] = (C[0] ^ C[2**N-1]) ^ (H[i][0] ^ H[i][1])
		G[i][1] = C[0] ^ H[i][0]
		#print i, C[0],  H[i][0],  H[i][1], (C[0] ^ C[2**N-1])
		#print "Res ", G[i][0], G[i][1]

	return G

def calcHilbertGenes(N):
	C = calC_r(N)
	H = calcH(C, N)
	G = calcG(C, H, N)

	return G

def printCcodeArray(array):
	string = ""

	if len(array.shape) == 2:
		string += "{"

		for num in array:
			tmpNum = 0x00000000000000000000000000000000
			for bit in num:
				tmpNum = tmpNum << 1
				tmpNum += bit

			string += "%i," % tmpNum

		string = string[:-1] + "}"

		return string

	elif len(array.shape) == 3:
		string += "{"

		for num in array:
			string += "{"

			for element in num:
				tmpNum = 0x00000000000000000000000000000000
				for bit in element:
					tmpNum = tmpNum << 1
					tmpNum += bit

				string += "%i," % tmpNum

			string = string[:-1] + "},"

		string = string[:-1] + "}"

		return string

	else:
		print "No... I'm not doing this..."

def printReverseCcodeArray(array):
	revArray = np.zeros(shape=(len(array)), dtype=np.uint64)
	string = ""

	if len(array.shape) == 2:
		i = 0
		for num in array:
			tmpNum = 0x00000000000000000000000000000000
			for bit in num:
				tmpNum = tmpNum << 1
				tmpNum += bit

			revArray[tmpNum] = i
			i += 1

		string += "{"

		for num in revArray:
			string += "%i," % num

		string = string[:-1] + "}"

		return string

	else:
		print "No... I'm not doing this..."


def main(argv):

	if len(argv) == 2:
		N = int(argv[1])
	else:
		print "Wrong arguments"
		return

	stringC = "{"
	stringRevC = "{"
	stringH = "{"

	for i in range(N):
		N = i + 1
		C = calC_r(N)
		G = calcHilbertGenes(N)

		print "uint32_t C%i[] = %s;" % (N, printCcodeArray(C))
		print "uint32_t revC%i[] = %s;" % (N, printReverseCcodeArray(C))
		print "uint32_t H%i[][2] = %s;" % (N, printCcodeArray(G))

		stringC += "(uint32_t*)&C%i," % N
		stringRevC += "(uint32_t*)&revC%i," % N
		stringH += "(uint32_t*)&H%i," % N

	stringC = stringC[:-1] + "}"
	stringRevC = stringRevC[:-1] + "}"
	stringH = stringH[:-1] + "}"

	print "uint32_t * C[] = %s;" % stringC
	print "uint32_t * revC[] = %s;" % stringRevC
	print "uint32_t * H[] = %s;" % stringH


if __name__ == "__main__":
	main(sys.argv)