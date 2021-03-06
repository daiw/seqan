// ==========================================================================
//                                index_sa_bpr
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

#ifndef CORE_TESTS_INDEX_SA_BPR_TEST_INDEX_SA_BPR_H_
#define CORE_TESTS_INDEX_SA_BPR_TEST_INDEX_SA_BPR_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/index.h>


using namespace seqan;

template<typename TText>
void testCodeD(TText text){

    typedef typename Value<TText>::Type TValue;
    typedef typename ValueSize<TValue>::Type TAlphabetSize;

    const TAlphabetSize ALPHABETSIZE = ValueSize<TValue>::VALUE + 1;

	for (short d = 1; d < 4; ++d) {
		long tmpModulo = pow(ALPHABETSIZE, (d - 1));
		int last = code_d(text, d, 0, ALPHABETSIZE);
		for (unsigned i = 1; i < length(text); ++i) {
			int current_rek = code_d(text, d, i, ALPHABETSIZE, tmpModulo, last);
			int current_it = code_d(text, d, i, ALPHABETSIZE);
			SEQAN_ASSERT_EQ(current_it, current_rek);
			last = current_it;
		}
	}
}

template<typename TText>
void compareSuffixArrays(TText text, unsigned short d) {

	typedef typename SAValue<TText>::Type TSA;
	String<TSA> sa;
	String<TSA> sa2;

	resize(sa, lengthSum(text));

#if (SEQAN_ENABLE_DEBUG || SEQAN_ENABLE_TESTING) && SEQAN_ENABLE_PARALLELISM
	const clock_t begin_time = omp_get_wtime();
	std::cout << "Computing BPR... ";
	std::cout.flush();
#endif

	createSuffixArray(sa, text, Bpr(), d);

#if (SEQAN_ENABLE_DEBUG || SEQAN_ENABLE_TESTING) && SEQAN_ENABLE_PARALLELISM
	std::cout << " done. " << float(omp_get_wtime() - begin_time);
	std::cout << std::endl;
	std::cout.flush();
#endif

	resize(sa2, lengthSum(text));

#if (SEQAN_ENABLE_DEBUG || SEQAN_ENABLE_TESTING) && SEQAN_ENABLE_PARALLELISM
	const clock_t begin_time2 = omp_get_wtime();
	std::cout << "Computing Skew3...   ";
	std::cout.flush();
#endif

	createSuffixArray(sa2, text, Skew3());

#if (SEQAN_ENABLE_DEBUG || SEQAN_ENABLE_TESTING) && SEQAN_ENABLE_PARALLELISM
	std::cout << " done. " << float(omp_get_wtime() - begin_time2);
	std::cout << std::endl;
	std::cout.flush();
#endif

	//Check for differences:
	unsigned errors = 0;
	for (unsigned i = 0; i < length(sa) && errors < 100; ++i) {
		SEQAN_ASSERT_EQ(sa[i], sa2[i]);
//		if (sa[i] != sa2[i]) {
//			std::cerr << "SuffixArray error at index " << i << std::endl;
//			++errors;
//		}
//		std::cerr.flush();
	}
//	if (errors == 0)
//		std::cout << "FEHLERFREI!" << std::endl;
}


// A test for code_D.
SEQAN_DEFINE_TEST(test_index_sa_bpr_codeD)
{

	using namespace seqan;

    String<Dna> text = "ACGTGCTG";
    testCodeD(text);

    String<char> text2 = "halloWelt!";
    testCodeD(text2);

    text = "AC";
    testCodeD(text);

    text = "";
    testCodeD(text);
}

// This test compares SuffixArrays Computed via Bpr with SuffixArrays computed via Skew
SEQAN_DEFINE_TEST(test_index_sa_bpr_compareSA)
{

	using namespace seqan;

    String<Dna> text = "ACGTGCTG";
    compareSuffixArrays(text, 5);

    String<char> text2 = "halloWelt!";
    compareSuffixArrays(text2, 3);

    text = "AC";
    compareSuffixArrays(text, 5);

    text = "";
    compareSuffixArrays(text, 5);
}

// This test compares SuffixArrays Computed via Bpr with SuffixArrays computed via Skew
SEQAN_DEFINE_TEST(test_index_sa_bpr_compareSAStringSets)
{
	using namespace seqan;

    StringSet<String<Dna> > text;

//    resize(text, 0);
//	compareSuffixArrays(text, 5);

    resize(text, 1);
    text[0] = "";
    compareSuffixArrays(text, 5);

    text[0] = "A";
    compareSuffixArrays(text, 5);

    text[0] = "ATGCAACGTVA";
    compareSuffixArrays(text, 5);

    resize(text, 2);
    text[1] = "";
    compareSuffixArrays(text, 5);

    text[1] = "A";
    compareSuffixArrays(text, 5);

    text[1] = "ATGCAACGTVA";
	compareSuffixArrays(text, 5);

	text[1] = "ACTAGCAGACGATAC";
	compareSuffixArrays(text, 5);

	StringSet<String<char> > text2;

	resize(text2, 1);
	text2[0] = "";
	compareSuffixArrays(text2, 3);

	text2[0] = "b";
	compareSuffixArrays(text2, 3);

	text2[0] = "halloWelt!";
	compareSuffixArrays(text2, 3);

	resize(text2, 2);
	text2[1] = "";
	compareSuffixArrays(text2, 3);

	text2[1] = "x";
	compareSuffixArrays(text2, 3);

	text2[1] = "halloWelt!";
	compareSuffixArrays(text2, 3);

	text2[1] = "halloWelt!12345";
	compareSuffixArrays(text2, 3);
}


#endif  // CORE_TESTS_INDEX_SA_BPR_TEST_INDEX_SA_BPR_H_
