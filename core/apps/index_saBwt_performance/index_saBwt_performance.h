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



#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>
#include <seqan/index.h>
#include <seqan/file.h>      // For printing SeqAn Strings.
#include <seqan/stream.h>

using namespace seqan;

// --------------------------------------------------------------------------
// Class AppOptions
// --------------------------------------------------------------------------
struct SaBwtAppOptions
{
    // Verbosity level.  0 -- quiet, 1 -- normal, 2 -- verbose, 3 -- very verbose.
    int verbosity;

    bool computeSA;

    bool computeBwt;

    bool compare;

    bool set;

    bool ascii;

    seqan::CharString file;

    SaBwtAppOptions() :
        verbosity(0), computeSA(false), computeBwt(false), compare(false), set(false), ascii(false)
    {}
};

template<typename TText, typename TSA, typename TTag>
void computeBprSuffixArray(TText& text, TSA& sa, const TTag tag) {

	resize(sa, lengthSum(text));
#ifdef _OPENMP
	const clock_t begin_time = omp_get_wtime();
#endif
	createSuffixArray(sa, text, tag);
#ifdef _OPENMP
	std::cout << " done in " << float(omp_get_wtime() - begin_time)<<"s. ";
#endif
}

template<typename TText, typename TSA>
void computeBprSuffixArray(TText& text, TSA& sa, Bpr const & tag) {

	resize(sa, lengthSum(text));

#ifdef _OPENMP
	const clock_t begin_time = omp_get_wtime();
#endif
	createSuffixArray(sa, text, tag, 5);
#ifdef _OPENMP
	std::cout << " done in " << float(omp_get_wtime() - begin_time)<<"s. ";
#endif
}

template<typename TSA>
void computeBprSuffixArray(CharString& text, TSA& sa, Bpr const & tag) {

    resize(sa, lengthSum(text));

#ifdef _OPENMP
    const clock_t begin_time = omp_get_wtime();
#endif
    createSuffixArray(sa, text, tag, 3);
#ifdef _OPENMP
    std::cout << " done in " << float(omp_get_wtime() - begin_time)<<"s. ";
#endif
}

template<typename TSA>
void compareSuffixArray(TSA& sa1, TSA& sa2) {

	//Check for differences:
	unsigned errors = 0;
	for (unsigned i = 0; i < length(sa1) && errors < 100; ++i) {
		if (sa1[i] != sa2[i]) {
			std::cerr << "SuffixArray error at index " << i << std::endl;
			++errors;
		}
		std::cerr.flush();
	}
	if (errors == 0)
		std::cout << "NO ERRORS!";

}

template<typename TText, typename TBWT, typename TSENT>
void computeBwt(StringSet<TText> &texts, TBWT &bwt, TSENT &sentinelPos) {

	typedef typename Value<TText>::Type TAlphabet;
	resize(bwt, lengthSum(texts) + countSequences(texts));
	resize(sentinelPos, countSequences(texts), Exact());
	TAlphabet sentinelChar = (TAlphabet) 0;

	float time = -1;
#ifdef _OPENMP
	const clock_t begin_time = omp_get_wtime();
#endif

	createBwt(bwt, sentinelPos, texts, sentinelChar);

#ifdef _OPENMP
	time = float(omp_get_wtime() - begin_time);
#endif
	std::cout << " done in " << time <<"s. ";
}

template<typename TText, typename TBWT, typename TSENT>
void computeBwt(TText &texts, TBWT &bwt, TSENT &sentinelPos) {

	typedef typename Value<TText>::Type TAlphabet;
	resize(bwt, lengthSum(texts) + countSequences(texts));
	resize(sentinelPos, countSequences(texts), Exact());
	TAlphabet sentinelChar = (TAlphabet) 0;

	float time = -1;
#ifdef _OPENMP
	const clock_t begin_time = omp_get_wtime();
#endif

	createBwt(bwt, sentinelPos, texts, sentinelChar);

#ifdef _OPENMP
	time = float(omp_get_wtime() - begin_time);
#endif
	std::cout << " done in " << time <<"s. ";
}

template<typename TText, typename TBWT, typename TSENT>
void computeBwtViaSA(StringSet<TText> &texts, TBWT &bwt, TSENT &sentinelPos) {

	typedef typename Value<TText>::Type TAlphabet;
	typedef typename SAValue<StringSet<TText> >::Type TSAValue;
	String<TSAValue> sa;
	resize(sa, lengthSum(texts), Exact());
	TAlphabet sentinelChar = (TAlphabet) 0;
	RankSupportBitString<void> dollarPos = _setDefaultSentinelPosition(
			length(bwt), RankSupportBitString<void>());
	resize(bwt, _computeBwtLength(texts), Exact());

	float time = -1;
#ifdef _OPENMP
	const clock_t begin_time = omp_get_wtime();
#endif

	createSuffixArray(sa, texts, Bpr(), 5);
	//createSuffixArray(sa, texts, Skew7());
	_createBwTable(bwt, dollarPos, texts, sa, sentinelChar);

#ifdef _OPENMP
	time = float(omp_get_wtime() - begin_time);
#endif
	std::cout << " done in " << time <<"s. ";

	resize(sentinelPos, countSequences(texts));

	int sentinelIndex = 0;
	for (unsigned i = 0; i < length(bwt); ++i) {
		if (bwt[i] == sentinelChar && isBitSet(dollarPos, i)) {
			sentinelPos[sentinelIndex++] = i;
		}
	}
}

template<typename TText, typename TBWT, typename TSENT>
void computeBwtViaSA(TText texts, TBWT &bwt, TSENT &sentinelPos) {

	typedef typename Value<TText>::Type TA;
	typedef typename Value<TA>::Type TAlphabet;
	typedef typename SAValue<TText>::Type TSAValue;
	String<TSAValue> sa;
	resize(sa, lengthSum(texts), Exact());

	TAlphabet sentinelChar = (TAlphabet) 0;
	unsigned dollarPos = 0;
	resize(bwt, _computeBwtLength(texts), Exact());

	float time = -1;
#ifdef _OPENMP
	const clock_t begin_time = omp_get_wtime();
#endif

	createSuffixArray(sa, texts, Bpr(), 5);
	//createSuffixArray(sa, texts, Skew7());
	_createBwTable(bwt, dollarPos, texts, sa, sentinelChar);

#ifdef _OPENMP
	time = float(omp_get_wtime() - begin_time);
#endif
	std::cout << " done in " << time <<"s. "<< std::endl;

	resize(sentinelPos, countSequences(texts));
	for (unsigned i = 0; i < length(bwt); ++i) {
		if (bwt[i] == sentinelChar && dollarPos==i) {
			sentinelPos[0] = i;
			break;
		}
	}
}

template<typename TBWT, typename TSENT>
void compareBwt(TBWT &bwt, TSENT &sentinelPos, TBWT &bwt2,
		TSENT &sentinelPos2) {

	if (length(bwt) != length(bwt2)) {
		std::cout << "Die beiden BWT sind nicht gleich lang: " << length(bwt)
				<< " vs. " << length(bwt2) << ".";
		std::cout << std::endl;
		return;
	}
	if (length(sentinelPos) != length(sentinelPos2)) {
		std::cout << "Die beiden SentinelLists sind nicht gleich lang."
				<< std::endl;
		return;
	}
	int error = 0;

	unsigned sentinelIndex = 0;
	for (unsigned i = 0; i < length(bwt); ++i) {
		if (sentinelIndex < length(sentinelPos)
				&& sentinelPos[sentinelIndex] == i) {
			if (sentinelPos[sentinelIndex] != sentinelPos2[sentinelIndex]) {
				std::cout << "Die Endzeichen stimmen nicht überein."
						<< std::endl;
				return;
			}
			++sentinelIndex;
		} else if (bwt[i] != bwt2[i]) {
			std::cout << "Zeichen an Stelle " << i << " ist falsch."
					<< std::endl;
			++error;
		}

		if(error>50)
			break;
	}
	if (error == 0)
		std::cout << "NO ERRORS! ";

}

template<typename TSA, typename TBwt, typename TInput>
void _internalDoTheWork(const SaBwtAppOptions& options, TInput& seqs) {

	int max = omp_get_max_threads();
	String<TSA> saRef;
	if (options.computeSA && options.compare) {
		std::cout << "\tComputing SA Skew7...";
		computeBprSuffixArray(seqs, saRef, Skew7());
		std::cout << std::endl;
	}
	TBwt bwtRef;
	String<unsigned> sentinelPosRef;
	if (options.computeBwt && options.compare) {
		std::cout << "\tComputing BWT via SA...";
		computeBwtViaSA(seqs, bwtRef, sentinelPosRef);
		std::cout << std::endl;
	}
	for (int k = max; k > 0; k--) {
		omp_set_num_threads(k);
		std::cout << "Processor count: " << omp_get_num_procs()
				<< " Thread count: " << omp_get_max_threads() << ". ";
		std::cout << std::endl;

		if (options.computeSA) {
			String<TSA> sa;
			std::cout << "\tComputing Bpr...";
			computeBprSuffixArray(seqs, sa, Bpr());
			if (options.compare) {
				compareSuffixArray(saRef, sa);
			}
			std::cout << std::endl;
		}
		if (options.computeBwt) {
			TBwt bwt;
			String<unsigned> sentinelPos;
			std::cout << "\tComputing Bwt...";
			computeBwt(seqs, bwt, sentinelPos);
			if (options.compare) {
				compareBwt(bwtRef, sentinelPosRef, bwt, sentinelPos);
			}
			std::cout << std::endl;
		}
	}
}

void doTheWork(SaBwtAppOptions options) {

    if(options.ascii){
        typedef String<char> TInput;
        typedef typename SAValue<TInput>::Type TSA;
        typedef String<char> TBwt;
        TInput buffer;

        std::fstream in(toCString(options.file), std::ios::in);
        typedef seqan::RecordReader<std::fstream, seqan::SinglePass<> > TRecordReader;
        TRecordReader reader(in);
       // readLetters(buffer, reader);

        while (!atEnd(reader))
        {
            char c = value(reader);
            appendValue(buffer, c, Generous());
            goNext(reader);
        }
        _internalDoTheWork<TSA, TBwt>(options, buffer);
    } else {
        typedef String<Dna5> TInput;
        typedef typename SAValue<TInput>::Type TSA;
        typedef String<Dna5> TBwt;
        TInput seqs;

        if (readFasta(seqs, toCString(options.file)) != 0){
            std::cout << "Error while reading fasta file";
            return;
        }
        _internalDoTheWork<TSA, TBwt>(options, seqs);
    }
}



void doTheWorkSet(SaBwtAppOptions options) {


	typedef StringSet<String<Dna5> > TInput;
	typedef typename SAValue<TInput>::Type TSA;
	typedef String<Dna5> TBwt;
	TInput seqs;

	std::fstream in(toCString(options.file),
		std::ios::binary | std::ios::in);
	seqan::RecordReader<std::fstream, seqan::SinglePass<> > reader(in);
	// Read file in one pass.
	{
		seqan::StringSet<seqan::CharString> ids;
		seqan::StringSet<seqan::CharString> quals;
		if (read2(ids, seqs, quals, reader, seqan::Fastq()) != 0){
		    std::cout << "Error while reading file";
			return;  // Could not read file.
		}
	}

	_internalDoTheWork<TSA, TBwt>(options, seqs);
}

