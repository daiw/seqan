// ==========================================================================
//                               index_bwt_bcr
// ==========================================================================
// Copyright (c) 2006-2013, Knut Reinert, FU Berlin
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Author: David Iwanowitsch <iwanowit@inf.fu-berlin.de>
// ==========================================================================

#ifndef CORE_TESTS_INDEX_BWT_BCR_TEST_INDEX_BWT_BCR_H_
#define CORE_TESTS_INDEX_BWT_BCR_TEST_INDEX_BWT_BCR_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/index.h>

// A test for strings.
SEQAN_DEFINE_TEST(test_index_bwt_bcr_sortBwtBucket)
{
    using namespace seqan;

	String<Triple<unsigned, char, bool> > bb;
	resize(bb, 8);
	String<Triple<unsigned, char, bool> > buffer;
	resize(buffer, 8);
	bb[0] = Triple<unsigned, char, bool>(0, 'a', false);
	bb[1] = Triple<unsigned, char, bool>(1, 'b', false);
	bb[2] = Triple<unsigned, char, bool>(2, 'f', false);
	bb[3] = Triple<unsigned, char, bool>(3, 'g', false);

	buffer[0] = Triple<unsigned, char, bool>(2, 'c', false);
	buffer[1] = Triple<unsigned, char, bool>(3, 'd', false);
	buffer[2] = Triple<unsigned, char, bool>(4, 'e', false);

	sortBwtBucket(bb, 4,  buffer, 3);

	for (unsigned i = 0; i < 7; ++i) {
		SEQAN_ASSERT_EQ(i, bb[i].i1);
		if(i>0){
			SEQAN_ASSERT_LT(bb[i - 1].i2, bb[i].i2);
		}
	}

	bb[0] =  Triple<unsigned, char, bool>(0, 'c', false);
	bb[1] =  Triple<unsigned, char, bool>(1, 'd', false);

	buffer[0] =  Triple<unsigned, char, bool>(0, 'a', false);
	buffer[1] =  Triple<unsigned, char, bool>(1, 'b', false);
	sortBwtBucket(bb, 2, buffer, 2);
	for (unsigned i = 0; i < 4; ++i) {
		SEQAN_ASSERT_EQ(i, bb[i].i1);
		if(i>0){
			SEQAN_ASSERT_LT(bb[i - 1].i2, bb[i].i2);
		}
	}

	bb[0] =  Triple<unsigned, char, bool>(0, 'a', false);
	bb[1] =  Triple<unsigned, char, bool>(1, 'b', false);

	buffer[0] =  Triple<unsigned, char, bool>(2, 'c', false);
	buffer[1] =  Triple<unsigned, char, bool>(3, 'd', false);
	sortBwtBucket(bb, 2, buffer, 2);
	for (unsigned i = 0; i < 4; ++i) {
		SEQAN_ASSERT_EQ(i, bb[i].i1);
		if(i>0){
			SEQAN_ASSERT_LT(bb[i - 1].i2, bb[i].i2);
		}
	}

	bb[0] =   Triple<unsigned, char, bool>(0, 'b', false);
	bb[1] =   Triple<unsigned, char, bool>(1, 'c', false);

	buffer[0] = Triple<unsigned, char, bool>(0, 'a', false);
	buffer[1] = Triple<unsigned, char, bool>(3, 'd', false);

	sortBwtBucket(bb, 2, buffer, 2);
	for (unsigned i = 0; i < 4; ++i) {
		SEQAN_ASSERT_EQ(i, bb[i].i1);
		if(i>0){
			SEQAN_ASSERT_LT(bb[i - 1].i2, bb[i].i2);
		}
	}

	bb[0] = Triple<unsigned, char, bool>(0, 'b', false);
	bb[1] = Triple<unsigned, char, bool>(1, 'c', false);
	bb[2] = Triple<unsigned, char, bool>(2, 'e', false);
	bb[3] = Triple<unsigned, char, bool>(3, 'f', false);

	buffer[0] = Triple<unsigned, char, bool>(0, 'a', false);
	buffer[1] = Triple<unsigned, char, bool>(3, 'd', false);

	sortBwtBucket(bb, 4, buffer, 2);
	for (unsigned i = 0; i < 6; ++i) {
		SEQAN_ASSERT_EQ(i, bb[i].i1);
		if(i>0){
			SEQAN_ASSERT_LT(bb[i - 1].i2, bb[i].i2);
		}
	}

	bb[0] = Triple<unsigned, char, bool>(0, 'b', false);
	bb[1] = Triple<unsigned, char, bool>(1, 'c', false);

	buffer[0] = Triple<unsigned, char, bool>(0, 'a', false);
	sortBwtBucket(bb, 2, buffer, 1);
	for (unsigned i = 0; i < 3; ++i) {
		SEQAN_ASSERT_EQ(i, bb[i].i1);
		if(i>0){
			SEQAN_ASSERT_LT(bb[i - 1].i2, bb[i].i2);
		}
	}
}

#endif  // CORE_TESTS_INDEX_BWT_BCR_TEST_INDEX_BWT_BCR_H_
