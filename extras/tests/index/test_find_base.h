// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2013 NVIDIA Corporation
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
//     * Neither the name of NVIDIA Corporation nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL NVIDIA CORPORATION BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Author: Enrico Siragusa <enrico.siragusa@fu-berlin.de>
// ==========================================================================

#ifndef EXTRAS_TESTS_FIND_BASE_H_
#define EXTRAS_TESTS_FIND_BASE_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/index.h>

using namespace seqan;

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Class Match
// ----------------------------------------------------------------------------

template <typename TTextOcc, typename TPatternOcc, typename TScore>
struct Match
{
    TTextOcc        textOcc;
    TPatternOcc     patternOcc;
    TScore          score;

    Match() :
        textOcc(0),
        patternOcc(0),
        score(0)
    {}

    Match(TTextOcc textOcc, TPatternOcc patternOcc, TScore score) :
        textOcc(textOcc),
        patternOcc(patternOcc),
        score(score)
    {}

    inline bool operator<(Match const & other) const
    {
        return (textOcc < other.textOcc) ||
        (textOcc == other.textOcc && patternOcc < other.patternOcc) ||
        (textOcc == other.textOcc && patternOcc == other.patternOcc && score < other.score);
    }

    inline bool operator>(Match const & other) const
    {
        return (textOcc > other.textOcc) ||
        (textOcc == other.textOcc && patternOcc > other.patternOcc) ||
        (textOcc == other.textOcc && patternOcc == other.patternOcc && score > other.score);
    }
};

// ----------------------------------------------------------------------------
// Class FinderTester
// ----------------------------------------------------------------------------

template <typename TText, typename TPattern, typename TSpec = void>
struct FinderTester
{
    typedef typename Fibre<TText, FibreSA>::Type                    TTextSAFibre;
    typedef typename Fibre<TPattern, FibreSA>::Type                 TPatternSAFibre;
    typedef typename Value<TTextSAFibre>::Type                      TTextSAPos;
    typedef typename Value<TPatternSAFibre>::Type                   TPatternSAPos;
    typedef Match<TTextSAPos, TPatternSAPos, unsigned>              TMatch;

    String<TMatch>   solution;
    String<TMatch>   results;
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function addSolution()                                        [FinderTester]
// ----------------------------------------------------------------------------

template <typename TText, typename TPattern, typename TSpec, typename TTextOcc, typename TPatternOcc, typename TScore>
inline void
addSolution(FinderTester<TText, TPattern, TSpec> & tester, TTextOcc textOcc, TPatternOcc patternOcc, TScore score)
{
    typedef FinderTester<TText, TPattern, TSpec>                                TTester;
    typedef typename TTester::TMatch                                            TMatch;

    appendValue(tester.solution, TMatch(textOcc, patternOcc, score));
}

// ----------------------------------------------------------------------------
// Function addResult()                                          [FinderTester]
// ----------------------------------------------------------------------------

template <typename TText, typename TPattern, typename TSpec, typename TTextOcc, typename TPatternOcc, typename TScore>
inline void
addResult(FinderTester<TText, TPattern, TSpec> & tester, TTextOcc textOcc, TPatternOcc patternOcc, TScore score)
{
    typedef FinderTester<TText, TPattern, TSpec>                                TTester;
    typedef typename TTester::TMatch                                            TMatch;

    appendValue(tester.results, TMatch(textOcc, patternOcc, score));
}

// ----------------------------------------------------------------------------
// Function test()                                               [FinderTester]
// ----------------------------------------------------------------------------

template <typename TText, typename TPattern, typename TSpec>
inline void
test(FinderTester<TText, TPattern, TSpec> & tester)
{
    std::sort(begin(tester.solution, Standard()), end(tester.solution, Standard()));
    std::sort(begin(tester.results, Standard()), end(tester.results, Standard()));

    SEQAN_ASSERT(isEqual(tester.solution, tester.results));
}

// ----------------------------------------------------------------------------
// Function buildTrie()
// ----------------------------------------------------------------------------
// TODO(esiragusa): move this into a Pattern class.

template <typename TIndex, typename TText>
void buildTrie(TIndex & index, TText const & text)
{
    typedef typename Iterator<TText const>::Type    TIterator;
    typedef typename Fibre<TIndex, FibreSA>::Type   TIndexSAFibre;
    typedef typename Value<TIndexSAFibre>::Type     TIndexSAPos;
    typedef typename Size<TText>::Type              TTextSize;
    
    TIndexSAFibre & sa = indexSA(index);
    reserve(sa, length(text), Exact());
    
    for (TIterator it = begin(text); it != end(text); it++)
    {
        TIndexSAPos localPos;
        assignValueI1(localPos, position(it));
        assignValueI2(localPos, 0);
        appendValue(sa, localPos, Exact());
    }

    QGramLess_<TIndexSAPos, TText const> less(text, MaxValue<TTextSize>::VALUE);
    std::sort(begin(sa, Standard()), end(sa, Standard()), less);
}

#endif  // EXTRAS_TESTS_FIND_BASE_H_
