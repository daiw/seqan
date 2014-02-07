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

#include <iostream>

#include <seqan/basic.h>
#include <seqan/index.h>
#include <seqan/file.h>      // For printing SeqAn Strings.

using namespace seqan;

int main(int argc, char const ** argv)
{
    StringSet<String<Dna> > texts;

    String<Dna> mtext1 = "TCTTA";
    appendValue(texts, mtext1);
    String<Dna> mtext2 = "CCCTA";
    appendValue(texts, mtext2);
    String<Dna> mtext3 = "CTTACTTA";
    appendValue(texts, mtext3);

    typedef typename SAValue<StringSet<String<Dna> > >::Type TSA;

    String<TSA> sa;
    resize(sa, lengthSum(texts));
    createSuffixArray(sa, texts, Bpr(), 5);

    std::cout << "BPR SA:   " << std::endl;
    for (unsigned i = 0; i < length(sa) && i < 800; ++i)
    {
        std::cout << sa[i] << " ";
        if (i > 0 && i % 100 == 0)
            std::cout << std::endl;
    }
    std::cout << std::endl;

    return 0;
}
