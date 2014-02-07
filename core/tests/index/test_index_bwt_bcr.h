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

using namespace seqan;

template<typename TText, typename TBWT, typename TSENT>
void computeBwt(TText &texts, TBWT &bwt, TSENT &sentinelPos)
{

    typedef typename Value<TText>::Type TA;
    typedef typename Value<TA>::Type TAlphabet;

    resize(bwt, lengthSum(texts) + countSequences(texts));
    resize(sentinelPos, countSequences(texts), Exact());
    TAlphabet sentinelChar = (TAlphabet) 0;

#if (SEQAN_ENABLE_DEBUG || SEQAN_ENABLE_TESTING) && SEQAN_ENABLE_PARALLELISM
    std::cout << "Computing BWT...   ";
    std::cout.flush();
    const clock_t begin_time = omp_get_wtime();
#endif

    createBwt(bwt, texts, sentinelPos, sentinelChar);

#if (SEQAN_ENABLE_DEBUG || SEQAN_ENABLE_TESTING) && SEQAN_ENABLE_PARALLELISM
    std::cout << " done. " << float(omp_get_wtime() - begin_time) << std::endl;
    std::cout.flush();
#endif
}

template<typename TText, typename TBWT, typename TSENT>
void computeBwtViaSA(StringSet<TText> &texts, TBWT &bwt, TSENT &sentinelPos)
{

    typedef typename Value<TText>::Type TAlphabet;
    typedef typename SAValue<StringSet<TText> >::Type TSAValue;
    String<TSAValue> sa;
    resize(sa, lengthSum(texts), Exact());
    TAlphabet sentinelChar = (TAlphabet) 0;
    RankSupportBitString<void> dollarPos = _setDefaultSentinelPosition(length(bwt), RankSupportBitString<void>());
    resize(bwt, _computeBwtLength(texts), Exact());

#ifdef _OPENMP
    std::cout << "Computing BWT via SA...   ";
    std::cout.flush();
    const clock_t begin_time = omp_get_wtime();
#endif

    createSuffixArray(sa, texts, Bpr(), 3);
    _createBwTable(bwt, dollarPos, texts, sa, sentinelChar);

#if (SEQAN_ENABLE_DEBUG || SEQAN_ENABLE_TESTING) && SEQAN_ENABLE_PARALLELISM
    std::cout << " done. " << float(omp_get_wtime() - begin_time) << std::endl;
    std::cout.flush();
#endif

    resize(sentinelPos, countSequences(texts));

    int sentinelIndex = 0;
    for (unsigned i = 0; i < length(bwt); ++i)
    {
        if (bwt[i] == sentinelChar && isBitSet(dollarPos, i))
        {
            sentinelPos[sentinelIndex++] = i;
        }
    }
}

template<typename TText, typename TBWT, typename TSENT>
void computeBwtViaSA(TText texts, TBWT &bwt, TSENT &sentinelPos)
{

    typedef typename Value<TText>::Type TA;
    typedef typename Value<TA>::Type TAlphabet;
    typedef typename SAValue<TText>::Type TSAValue;
    String<TSAValue> sa;
    resize(sa, lengthSum(texts), Exact());

    TAlphabet sentinelChar = (TAlphabet) 0;
    unsigned dollarPos = 0;
    resize(bwt, _computeBwtLength(texts), Exact());

#if (SEQAN_ENABLE_DEBUG || SEQAN_ENABLE_TESTING) && SEQAN_ENABLE_PARALLELISM
    std::cout << "Computing BWT via SA...  ";
    std::cout.flush();
    const clock_t begin_time = omp_get_wtime();
#endif

    createSuffixArray(sa, texts, Bpr(), 3);
    _createBwTable(bwt, dollarPos, texts, sa, sentinelChar);

#if (SEQAN_ENABLE_DEBUG || SEQAN_ENABLE_TESTING) && SEQAN_ENABLE_PARALLELISM
    std::cout << " done. " << float(omp_get_wtime() - begin_time) << std::endl;
#endif

    resize(sentinelPos, countSequences(texts));
    for (unsigned i = 0; i < length(bwt); ++i)
    {
        if (bwt[i] == sentinelChar && dollarPos == i)
        {
            sentinelPos[0] = i;
            break;
        }
    }
}

template<typename TBWT, typename TSENT>
void compareBwt(TBWT &bwt, TSENT &sentinelPos, TBWT &bwt2, TSENT &sentinelPos2)
{

    SEQAN_ASSERT_EQ_MSG(length(bwt), length(bwt2), "The bwts have a different length.");

    SEQAN_ASSERT_EQ_MSG(length(sentinelPos), length(sentinelPos2), "The sentinel lists have a different length.");

    unsigned sentinelIndex = 0;
    for (unsigned i = 0; i < length(bwt); ++i)
    {
        if (sentinelIndex < length(sentinelPos) && sentinelPos[sentinelIndex] == i)
        {

            SEQAN_ASSERT_EQ_MSG(sentinelPos[sentinelIndex], sentinelPos2[sentinelIndex],
                    "The bwts have different sentinel positions.");
            ++sentinelIndex;
        }
        else
        {
            SEQAN_ASSERT_EQ(bwt[i], bwt2[i]); //, "The bwts are different at position %d.", i);
        }
    }
}

// A test for comprating bwt_bcr with a bwt computed via SuffixArray.
SEQAN_DEFINE_TEST(test_index_bwt_bcr_compareBwt){
    using namespace seqan;

    {
        String<Dna> text1 = "TGCCAAC";
        String<Dna> text2 = "AGAGCTC";
        String<Dna> text3 = "GTCGCTT";

        StringSet<String<Dna> > texts;
        appendValue(texts, text1); appendValue(texts, text2); appendValue(texts, text3);

        String<Dna> bwt;
        String<unsigned> sentinelPos;
        String<Dna> bwt2;
        String<unsigned> sentinelPos2;
        computeBwt(texts, bwt, sentinelPos);
        computeBwtViaSA(texts, bwt2, sentinelPos2);
        compareBwt(bwt, sentinelPos, bwt2, sentinelPos2);
    }

    {
        String<char> bwtChar;
        String<unsigned> sentinelPosChar;
        String<char> bwtChar2;
        String<unsigned> sentinelPosChar2;
        String<char> text = "apple";

        computeBwt(text, bwtChar, sentinelPosChar);
        computeBwtViaSA(text, bwtChar2, sentinelPosChar2);
        compareBwt(bwtChar, sentinelPosChar, bwtChar2, sentinelPosChar2);
    }

    {
        String<char> bwtCharX;
        String<unsigned> sentinelPosCharX;
        String<char> bwtCharX2;
        String<unsigned> sentinelPosCharX2;
        String<char> textX2 = "BANANA";
        String<char> textX3 = "BANANAAsakjdgbasfislufkjds";
        String<char> textX4 = "BANANAfalsukdjfbdsavdxv";
        String<char> textX5 = "falsukjdgbfdasBANANA";
        String<char> textX6 = "asdfBANafsdlfkhjsfasANA";
        String<char> textX7 = "qa,jshfgs"; //
        StringSet<String<char> > textsX;
        appendValue(textsX, textX2);
        appendValue(textsX, textX3);
        appendValue(textsX, textX4);
        appendValue(textsX, textX5);
        appendValue(textsX, textX6);
        appendValue(textsX, textX7);
        computeBwt(textsX, bwtCharX, sentinelPosCharX);
        computeBwtViaSA(textsX, bwtCharX2, sentinelPosCharX2);
        compareBwt(bwtCharX, sentinelPosCharX, bwtCharX2, sentinelPosCharX2);
    }

    {
        String<Dna> text1 = "";
        String<Dna> text2 = "";
        String<Dna> text3 = "";

        StringSet<String<Dna> > texts;
        appendValue(texts, text1); appendValue(texts, text2); appendValue(texts, text3);

        String<Dna> bwt;
        String<unsigned> sentinelPos;
        String<Dna> bwt2;
        String<unsigned> sentinelPos2;
        computeBwt(texts, bwt, sentinelPos);
        //bwt via SA crashes with empty strings
    //		computeBwtViaSA(texts, bwt2, sentinelPos2);
    //		compareBwt(bwt, sentinelPos, bwt2, sentinelPos2);
    }

    {
        String<Dna> text1 = "";
        String<Dna> text2 = "AGAGCTC";
        String<Dna> text3 = "";

        StringSet<String<Dna> > texts;
        appendValue(texts, text1); appendValue(texts, text2); appendValue(texts, text3);

        String<Dna> bwt;
        String<unsigned> sentinelPos;
        String<Dna> bwt2;
        String<unsigned> sentinelPos2;
        computeBwt(texts, bwt, sentinelPos);
        //bwt via SA crashes with empty strings
        //computeBwtViaSA(texts, bwt2, sentinelPos2);
        //compareBwt(bwt, sentinelPos, bwt2, sentinelPos2);
    }

    {
        String<Dna> text1 = "";
        String<Dna> text2 = "";
        String<Dna> text3 = "GTCGCTT";

        StringSet<String<Dna> > texts;
        appendValue(texts, text1); appendValue(texts, text2); appendValue(texts, text3);

        String<Dna> bwt;
        String<unsigned> sentinelPos;
        String<Dna> bwt2;
        String<unsigned> sentinelPos2;
        computeBwt(texts, bwt, sentinelPos);
        //bwt via SA crashes with empty strings
    //		computeBwtViaSA(texts, bwt2, sentinelPos2);
    //		compareBwt(bwt, sentinelPos, bwt2, sentinelPos2);
    }

    {
        String<Dna> text1 = "TGCCAAC";
        String<Dna> text2 = "";
        String<Dna> text3 = "";

        StringSet<String<Dna> > texts;
        appendValue(texts, text1); appendValue(texts, text2); appendValue(texts, text3);

        String<Dna> bwt;
        String<unsigned> sentinelPos;
        String<Dna> bwt2;
        String<unsigned> sentinelPos2;
        computeBwt(texts, bwt, sentinelPos);
        //bwt via SA crashes with empty strings
    //		computeBwtViaSA(texts, bwt2, sentinelPos2);
    //		compareBwt(bwt, sentinelPos, bwt2, sentinelPos2);
    }

    {
        String<Dna> text1 = "T";
        String<Dna> text2 = "AGAGCTC";
        String<Dna> text3 = "";

        StringSet<String<Dna> > texts;
        appendValue(texts, text1); appendValue(texts, text2); appendValue(texts, text3);

        String<Dna> bwt;
        String<unsigned> sentinelPos;
        String<Dna> bwt2;
        String<unsigned> sentinelPos2;
        computeBwt(texts, bwt, sentinelPos);
        //bwt via SA crashes with empty strings
    //		computeBwtViaSA(texts, bwt2, sentinelPos2);
    //		compareBwt(bwt, sentinelPos, bwt2, sentinelPos2);
    }

}

// A test for sortBwtBucket.
SEQAN_DEFINE_TEST(test_index_bwt_bcr_sortBwtBucket){
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

    sortBwtBucket(bb, 4, buffer, 3);

    for (unsigned i = 0; i < 7; ++i)
    {
        SEQAN_ASSERT_EQ(i, bb[i].i1);
        if(i>0)
        {
            SEQAN_ASSERT_LT(bb[i - 1].i2, bb[i].i2);
        }
    }

    bb[0] = Triple<unsigned, char, bool>(0, 'c', false);
    bb[1] = Triple<unsigned, char, bool>(1, 'd', false);

    buffer[0] = Triple<unsigned, char, bool>(0, 'a', false);
    buffer[1] = Triple<unsigned, char, bool>(1, 'b', false);
    sortBwtBucket(bb, 2, buffer, 2);
    for (unsigned i = 0; i < 4; ++i)
    {
        SEQAN_ASSERT_EQ(i, bb[i].i1);
        if(i>0)
        {
            SEQAN_ASSERT_LT(bb[i - 1].i2, bb[i].i2);
        }
    }

    bb[0] = Triple<unsigned, char, bool>(0, 'a', false);
    bb[1] = Triple<unsigned, char, bool>(1, 'b', false);

    buffer[0] = Triple<unsigned, char, bool>(2, 'c', false);
    buffer[1] = Triple<unsigned, char, bool>(3, 'd', false);
    sortBwtBucket(bb, 2, buffer, 2);
    for (unsigned i = 0; i < 4; ++i)
    {
        SEQAN_ASSERT_EQ(i, bb[i].i1);
        if(i>0)
        {
            SEQAN_ASSERT_LT(bb[i - 1].i2, bb[i].i2);
        }
    }

    bb[0] = Triple<unsigned, char, bool>(0, 'b', false);
    bb[1] = Triple<unsigned, char, bool>(1, 'c', false);

    buffer[0] = Triple<unsigned, char, bool>(0, 'a', false);
    buffer[1] = Triple<unsigned, char, bool>(3, 'd', false);

    sortBwtBucket(bb, 2, buffer, 2);
    for (unsigned i = 0; i < 4; ++i)
    {
        SEQAN_ASSERT_EQ(i, bb[i].i1);
        if(i>0)
        {
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
    for (unsigned i = 0; i < 6; ++i)
    {
        SEQAN_ASSERT_EQ(i, bb[i].i1);
        if(i>0)
        {
            SEQAN_ASSERT_LT(bb[i - 1].i2, bb[i].i2);
        }
    }

    bb[0] = Triple<unsigned, char, bool>(0, 'b', false);
    bb[1] = Triple<unsigned, char, bool>(1, 'c', false);

    buffer[0] = Triple<unsigned, char, bool>(0, 'a', false);
    sortBwtBucket(bb, 2, buffer, 1);
    for (unsigned i = 0; i < 3; ++i)
    {
        SEQAN_ASSERT_EQ(i, bb[i].i1);
        if(i>0)
        {
            SEQAN_ASSERT_LT(bb[i - 1].i2, bb[i].i2);
        }
    }
}

#endif  // CORE_TESTS_INDEX_BWT_BCR_TEST_INDEX_BWT_BCR_H_
