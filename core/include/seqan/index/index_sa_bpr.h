// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
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

#ifndef CORE_INCLUDE_SEQAN_INDEX_INDEX_SA_BPR_H_
#define CORE_INCLUDE_SEQAN_INDEX_INDEX_SA_BPR_H_

namespace seqan
{

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

struct Bpr_;
typedef Tag<Bpr_> Bpr;

// comparator used for storting
// for a given index it returns the BucketPointerValue for the suffixArray entry,
// while respecting the offsets within StringSets
template<typename InType, typename TSa, typename TBptr, typename TLimits, typename TOffset,
        typename OutType = typename Value<TBptr>::Type>
struct _bprComparator : public std::unary_function<InType, OutType>
{
    const TSa &SA;
    const TBptr &bptr;
    const TLimits &limits;
    const long bptrExtPerString;
    const TOffset offset;

    _bprComparator(const TSa &_sa, const TBptr &_bptr, const TLimits &_limits, const long _bptrExtPerString,
            const TOffset _offset) :
            SA(_sa), bptr(_bptr), limits(_limits), bptrExtPerString(_bptrExtPerString), offset(_offset)
    {
    }

    OutType operator()(InType index)
    {
        return getBptrVal(SA, bptr, limits, bptrExtPerString, offset, index);
    }
};

// ============================================================================
// Metafunctions
// ============================================================================

/**
.Function.createSuffixArray:
..summary:Computes the SuffixArray for given input
..description:
Creates a suffixarray for the input text. The algorithm is based on the bucket pointer refinement algorithm by
Schürmann and Stoye. It is modified and extended for parallel computation.
..createSuffixArray(SA, text, tag, qgram-length)
..param.qgram-length: Suffixlength for initial bucket sort.
*/
template<typename TSA, typename TText>
void createSuffixArray(TSA & SA, TText & s, Bpr const &, unsigned short const d)
{
    typedef typename Value<TSA>::Type TSAValue;
    typedef typename Size<TText>::Type TTextSize;
    typedef typename StringSetLimits<TText>::Type TStringSetLimits;

    if (lengthSum(s) < 1)
        return;

    const TStringSetLimits limits = stringSetLimits(s);

/////////////////////////
// Phase 1:
// sort suffixes with length d and create bucketPointer bptr
//
/////////////////////////

    String<TTextSize> bkt; //set contains the count of the alphabetsize^d buckets

    //contains for each entry in text, the index of the last entry of the corresponding Bucket in SA
    //remark: contains negative values to get the correct order of the $-signs
    String<long> bptr;

#if (SEQAN_ENABLE_DEBUG || SEQAN_ENABLE_TESTING) && SEQAN_ENABLE_PARALLELISM
    const double begin_time = omp_get_wtime();
#endif

    fillBucketsPhase1(SA, bptr, bkt, s, limits, d, getAlphabetSize(s)+1);

#if (SEQAN_ENABLE_DEBUG || SEQAN_ENABLE_TESTING) && SEQAN_ENABLE_PARALLELISM
    std::cout << "\t Phase 1: " << double(omp_get_wtime() - begin_time)
    << "s. ";
#endif

/////////////////////////
// Phase 2:
// recursively refine buckets with size > 1 by sorting buckets with corresponding
// offset and refresh bptr
/////////////////////////

#if (SEQAN_ENABLE_DEBUG || SEQAN_ENABLE_TESTING) && SEQAN_ENABLE_PARALLELISM
    const double begin_time_2 = omp_get_wtime();
#endif

    sewardCopyPhase2(SA, bptr, s, bkt, limits, d, getAlphabetSize(s)+1);

    clear(bkt);
    clear(bptr);

#if (SEQAN_ENABLE_DEBUG || SEQAN_ENABLE_TESTING) && SEQAN_ENABLE_PARALLELISM
    std::cout << "\t Phase 2: "
    << double(omp_get_wtime() - begin_time_2) << "s. ";
#endif

}

/**
.Function.getAlphabetSize:
..summary:returns the Alphabetsize of the String
*/
template<typename TString>
unsigned getAlphabetSize(TString const &){
    typedef typename Value<TString>::Type ARGH;
    return ValueSize<ARGH>::VALUE;
}

/**
.Function.getAlphabetSize:
..summary:returns the Alphabetsize of the inner String
*/
template<typename TString>
unsigned getAlphabetSize(StringSet<TString> const &){
    typedef typename Value<TString>::Type ARGH;
    return ValueSize<ARGH>::VALUE;
}

/*
 * sort suffixes with length d and create Bucketpointer
 *
 * Set bucketpointer to buckets last position in SA
 */
template<typename TSA, typename TText, typename TBptr, typename TBkt, typename TLimit, typename TD, typename TAlpha>
inline void fillBucketsPhase1(TSA & SA, TBptr & bptr, TBkt & bkt, TText const & s, const TLimit & limits, const TD d,
        const TAlpha alphabeSizeWithDollar)
{
    typedef typename SAValue<TText>::Type TSAValue;
    typedef typename LengthSum<TText>::Type TTextSize;
    typedef typename Size<TText>::Type TTextCount;
    typedef typename Value<TBkt>::Type TBktValue;
    typedef typename Size<TBkt>::Type TBktSize;
    typedef typename Size<Splitter<TTextSize> >::Type TSplitterSize;

    unsigned const bptrExtPerString = (2 * d + 1);//after each sequence we insert 2*d+1 values to sort the $ correctly

    //precompute some values:
    unsigned const tmpModulo = pow(alphabeSizeWithDollar, (d - 1));
    unsigned const bucketCount = tmpModulo * alphabeSizeWithDollar;

    TTextSize const n = lengthSum(s);
    TTextCount const stringCount = countSequences(s);

    //bptr has the size of n plus an offset for each string in a stringSet
    resize(bptr, n + stringCount * bptrExtPerString, Exact());
    resize(bkt, bucketCount + 1, Exact());

    //Splitt inputlength to equal parts per thread
    Splitter<TTextSize> splitter(0, n);

    //each thread gets his own
    String<TBkt> bktPerThread;
    resize(bktPerThread, length(splitter), Exact());

    //parallel: Split text in parts of equal size for each thread
    SEQAN_OMP_PRAGMA(parallel num_threads(length(splitter)))
    {
        unsigned const threadNum = omp_get_thread_num();
        TTextSize const start = splitter[threadNum];
        TTextSize const end = splitter[threadNum + 1];
        TSAValue saValue;
        long bptrOffset;

        //init Bucket per Thread
        TBkt &threadBucket = bktPerThread[threadNum];
        resize(threadBucket, bucketCount + 1, 0, Exact());

        //first iteration: compute rank and bucketsizes
        posLocalize(saValue, start, limits);
        bptrOffset = getSeqNo(saValue) * bptrExtPerString;

        unsigned codeD = code_d(
                getSequenceByNo(getSeqNo(saValue), s), d, saValue,
                alphabeSizeWithDollar);
        bptr[start + bptrOffset] = codeD;
        threadBucket[codeD]++;

        for (TTextSize i = start + 1; i < end; ++i)
        {
            posLocalize(saValue, i, limits);
            bptrOffset = getSeqNo(saValue) * bptrExtPerString;

            codeD = code_d(
                    getSequenceByNo(getSeqNo(saValue), s), d, saValue,
                    alphabeSizeWithDollar, tmpModulo,
                    codeD);
            bptr[i + bptrOffset] = codeD;

            threadBucket[codeD]++;
        }

        SEQAN_OMP_PRAGMA(barrier)

        //summarize buckets: each thread needs to know how many entries the previous threads are going to insert into the bucket
        //bkt contains the starting positions (=number of entrys in smaller buckets)
        SEQAN_OMP_PRAGMA(single)
        {
            bkt[0] = 0;
            TBktValue last = 0;
            TBktValue tmp;
            for (TBktSize i = 0; i < length(bktPerThread[0]) - 1; ++i)
            {
                if (i == bucketCount)
                    continue;

                //we fill the next bucketentry with the previous count, to have the bucket-startposition
                bkt[i + 1] = last;

                //We do the same per thread, so that each thread has its own startposition and we avoid locking
                TBktValue perThreadLast = 0;
                for (TSplitterSize j = 0; j < length(splitter); ++j)
                {
                    tmp = bktPerThread[j][i];
                    bktPerThread[j][i] = perThreadLast;
                    perThreadLast += tmp;
                }
                bkt[i + 1] += perThreadLast;

                last = bkt[i + 1];
            }
            bkt[bucketCount] = n;
        }

        //now fill the buckets
        for (TTextSize i = start; i < end; ++i)
        {
            posLocalize(saValue, i, limits);
            bptrOffset = getSeqNo(saValue) * bptrExtPerString;
            TBktValue index = (bkt[bptr[i + bptrOffset] + 1]
                    - threadBucket[bptr[i + bptrOffset]]) - 1;

            threadBucket[bptr[i + bptrOffset]]++;
            SA[index] = saValue;
        }

        //up to now we misused bptr to save the rank-value of a suffix
        //now bptr should be a pointer from suffix to bucket
        for (TTextSize i = start; i < end; ++i)
        {
            posLocalize(saValue, i, limits);
            bptrOffset = getSeqNo(saValue) * bptrExtPerString;
            bptr[i + bptrOffset] = bkt[bptr[i + bptrOffset] + 1] - 1;
        }

        SEQAN_OMP_PRAGMA(barrier)

        //the last d-suffixes get special handling corresponding to the sorting of the $-signs
        SEQAN_OMP_PRAGMA(single)
        {
            for (TTextCount k = 1; k <= stringCount; ++k)
            {
                TTextCount seq = stringCount - k;
                bptrOffset = seq * bptrExtPerString;

                TTextSize const seqEndIndex = sequenceLength(seq, s);
                TTextSize const seqStartIndex = posGlobalize(
                        Pair<TTextCount, TTextSize>(seq, 0), limits);
                unsigned const tmpValue = code_d(getSequenceByNo(seq, s), d,
                        seqEndIndex - d - 1, alphabeSizeWithDollar);
                unsigned lastValue = tmpValue;

                for (TTextSize i = seqEndIndex - d; i < seqEndIndex; ++i)
                {
                    lastValue = code_d(getSequenceByNo(seq, s), d, i,
                            alphabeSizeWithDollar, tmpModulo, lastValue);

                    bptr[bptrOffset + seqStartIndex + i] = bkt[lastValue];
                    bkt[lastValue] = bkt[lastValue] + 1;
                }

                //insert negative values after each sequence
                TTextSize j = seq * 2 * d + seq;
                j *= -1;
                for (TTextSize i = seqEndIndex; i <= seqEndIndex + 2 * d; ++i)
                {
                    bptr[seqStartIndex + i + bptrOffset] = --j;
                }
            }

            //reinsert correct bkt values (needed for seward copy)
            for (TTextCount k = 1; k <= stringCount; ++k)
            {
                TTextCount seq = stringCount - k;
                bptrOffset = seq * bptrExtPerString;

                TTextSize const seqEndIndex = sequenceLength(seq, s);
                unsigned lastValue = code_d(getSequenceByNo(seq, s), d,
                        seqEndIndex - d - 1, alphabeSizeWithDollar);
                for (TTextSize i = seqEndIndex - d; i < seqEndIndex; ++i)
                {
                    lastValue = code_d(getSequenceByNo(seq, s), d, i,
                            alphabeSizeWithDollar, tmpModulo, lastValue);
                    bkt[lastValue] = bkt[lastValue] - 1;
                }
            }
        }
    }

    clear(bktPerThread);
}

// computes the "multiple character encoding" (the rank)
//
//&s: SuffixArray
//d: length
//i: pos in s
template<typename TText, typename TPos, typename TAlphabetSize>
unsigned code_d(TText const & s, unsigned short const & d, TPos const & pos, TAlphabetSize const & alphabetSize)
{
    unsigned int result = 0;
    unsigned local = getSeqOffset(pos);

    for (unsigned short k = 1; k <= d; ++k)
    {
        if (length(s) <= (local + k - 1))
            result += 0;
        else
            result += pow(alphabetSize, (d - k)) * (unsigned) (1 + ordValue(s[local + k - 1]));
    }
    return result;
}

// computes the "multiple character encoding" using the previous result
//
//&s: SuffixArray
//d: length
//i: pos in s
//code_d_i: value of previous suffix
template<typename TText, typename TPos, typename TAlphabetSize>
unsigned code_d(TText const & s, unsigned short const & d, TPos const & pos, TAlphabetSize const & alphabetSize, unsigned const & modulo,
        unsigned const & code_d_i)
{
    unsigned local = getSeqOffset(pos);
    if (local == 0)
        return code_d(s, d, pos, alphabetSize);

    if (length(s) <= local + d - 1)
        return alphabetSize * (code_d_i % modulo) + 0;
    else
        return alphabetSize * (code_d_i % modulo) + 1 + ordValue(s[local + d - 1]);
}

//Phase 2: sort and refine all Buckets
template<typename TSA, typename TText, typename TBptr, typename TBkt, typename TLimit, typename TD, typename TAlpha>
inline void sewardCopyPhase2(TSA & SA, TBptr & bptr, TText const & s, TBkt const & bkt, TLimit const & limits, TD const d,
        TAlpha const alphabeSizeWithDollar)
{
    typedef typename Value<TSA>::Type TSAValue;
    typedef typename Size<TText>::Type TTextSize;
    typedef typename Value<TBptr>::Type TBptrValue;
    typedef typename Value<TBkt>::Type TBktIndex;

    const TAlpha ALPHABETSIZE = getAlphabetSize(s);
    const unsigned bptrExtPerString = (2 * d + 1);

    const unsigned bucketsInL2Bucket = pow(alphabeSizeWithDollar, d - 2); //d minus 2, since these chars are equal in L2
    const unsigned bucketsInL1Bucket = bucketsInL2Bucket * alphabeSizeWithDollar;

    //A String for each Thread containing Triple(startIndex, endIndex, offset)
    StringSet<String<Triple<TBktIndex, TBktIndex, unsigned> > > bucketIndicesSet;
    resize(bucketIndicesSet, omp_get_max_threads());

    //Initial: find non empty Buckets and store their indices
    SEQAN_OMP_PRAGMA(parallel)
    {
        String<Triple<TBktIndex, TBktIndex, unsigned> > &bucketIndicesCur = bucketIndicesSet[omp_get_thread_num()];
        SEQAN_OMP_PRAGMA(for collapse(2) schedule(dynamic))
        for (TAlpha i = 0; i < ALPHABETSIZE; ++i)
        {
            for (TAlpha j = 0; j < ALPHABETSIZE; ++j)
            {
                //only take buckets, starting with two different characters to avoid inefficient refinement of large 'NNN...N' sequences
                if (j == i)
                    continue;

                //add 1 since we use 0 as the dollar-sign
                const TAlpha firstChar = 1 + i;
                const TAlpha secondChar = 1 + j;

                const long bktStartIndex = firstChar * bucketsInL1Bucket + secondChar * bucketsInL2Bucket; //startbucket for suffixes starting with firstchar, secondchar
                const long bktEndIndex = bktStartIndex + bucketsInL2Bucket;

                for (long k = bktStartIndex; k < bktEndIndex; ++k)
                {
                    //all buckets starting with firstchar, secondchar
                    const TBktIndex left = bkt[k]; //Start
                    const TBktIndex right = bkt[k + 1] - 1; //end: next bucket start -1

                    long bucketSize = right - left;
                    if (bucketSize > 0)
                    {
                        //append this Bucket to list
                        appendValue(bucketIndicesCur, Triple<TBktIndex, TBktIndex, unsigned>(left, right, d));
                    }
                }
            }
        }
    }

#if (SEQAN_ENABLE_DEBUG || SEQAN_ENABLE_TESTING) && SEQAN_ENABLE_PARALLELISM
    const clock_t begin_time = omp_get_wtime();
#endif

    /*
     * Now sort these Buckets
     */

    StringSet<String<Triple<TBktIndex, TBktIndex, unsigned> > > nextBucketIndices;
    String<Triple<TBktIndex, TBktIndex, unsigned> > bucketIndices = concat(bucketIndicesSet);
    String<long> split;
    String<Pair<TBktIndex, TBktIndex> > resultValues;

    String<String<Pair<TBktIndex, TBktIndex> > > nextBptrValuesPerThread;
    resize(nextBptrValuesPerThread, omp_get_max_threads());

    //minimal bucketcount of each job
    const long minBuckets = 100 * omp_get_max_threads();

    while (length(bucketIndices) > 0)
    {
        //Split buckets in reasonable chunks, after each chunk sync refined bptr of each thread
        clear(split);
        appendValue(split, 0);
        long currentLength = 0;
        long lastSplit = 0;
        for (unsigned k = 0; k < length(bucketIndices); ++k)
        {
            if (k - lastSplit > minBuckets)
            {
                appendValue(split, k);
                currentLength = 0;
                lastSplit = k;
            }
            currentLength += bucketIndices[k].i2 - bucketIndices[k].i1;
        }
        appendValue(split, length(bucketIndices));
        resize(resultValues, length(bucketIndices));

        //sort each split and refine afterwards in parallel
        //we can not sort and refine in parallel
        //if we wait to long with the refinement, the sorting get uneffective.
        //so we use split to regulary switch between sorting and refinement.
        for (unsigned i = 1; i < length(split); ++i)
        {
            const long startIndex = split[i - 1];
            const long endIndex = split[i];

            SEQAN_OMP_PRAGMA(parallel shared(bucketIndices, nextBucketIndices, SA, bptr, limits))
            {
                //sort all collected buckets in parallel
                SEQAN_OMP_PRAGMA(for schedule(guided))
                for (unsigned k = startIndex; k < endIndex; ++k)
                {
                    const unsigned offset = bucketIndices[k].i3;
                    const TBktIndex start = bucketIndices[k].i1;
                    const TBktIndex end = bucketIndices[k].i2;

                    if(end - start == 1)
                    {
                        sortSizeTwo(SA, bptr, limits, d, offset, start, end);
                        bucketIndices[k].i3 = 0; //marker to ignore this bucket after refinement (no recursion needed)
                    }
                    else
                    {
                        typedef _bprComparator<TTextSize, TSA, TBptr, TLimit, TTextSize> TSortFunctor;
                        TSortFunctor sortFunctor = TSortFunctor(SA, bptr, limits, bptrExtPerString, offset);
                        doQuickSort(QsortSequential(), SA, sortFunctor, start, end);
                    }
                }

                //refine the previously sorted buckets in parallel

                //collect all changed bptr values per thread:
                String<Pair<TBktIndex, TBktIndex> > &nextBptrValues = nextBptrValuesPerThread[omp_get_thread_num()];

                SEQAN_OMP_PRAGMA(for schedule(guided))
                for (unsigned k = startIndex; k < endIndex; ++k)
                {
                    const unsigned offset = bucketIndices[k].i3;
                    const TBktIndex start = bucketIndices[k].i1;
                    const TBktIndex end = bucketIndices[k].i2;

                    //TODO: erkläre resultValues
                    resultValues[k] = refineBucket(SA, bptr, nextBptrValues, limits, start, end, offset, d);
                }

                // note, the implicit barrier after the for-loop above

                //sync bptr, to use refined pointers in next chunk
                for (unsigned k = 0; k < length(nextBptrValues); ++k)
                {
                    bptr[nextBptrValues[k].i1]= nextBptrValues[k].i2;
                }

                clear(nextBptrValues);
            }
        }

        resize(nextBucketIndices, omp_get_max_threads());

        //search next buckets for next iteration in parallel
        SEQAN_OMP_PRAGMA(parallel for schedule(guided))
        for (unsigned k = 0; k < length(bucketIndices); ++k)
        {
            const unsigned offset = bucketIndices[k].i3;
            if(offset == 0)
                continue;
            const TBktIndex start = bucketIndices[k].i1;
            const TBktIndex end = bucketIndices[k].i2;

            unsigned bptrOffset = getSeqNo(SA[start]) * bptrExtPerString;

            const Pair<TBktIndex, TBktIndex> middleValues = resultValues[k];

            String<Triple<TBktIndex, TBktIndex, unsigned> > &nextBucketIndicesThread = nextBucketIndices[omp_get_thread_num()];

            //increase offset
            unsigned int newOffset = offset + d;
            if ((TBktIndex)bptr[bptrOffset + posGlobalize(SA[start], limits)] == end)
            {
                //further increase offset
                newOffset = computeLCPAndOffset(SA, bptr, limits, start, end, newOffset, d);
            }

            //TODO: das noch ein bisschen kommentieren
            TBktIndex leftTmp = start;
            while (leftTmp < middleValues.i1)
            {
                bptrOffset = getSeqNo(SA[leftTmp]) * bptrExtPerString;
                const TBktIndex rightTmp = bptr[bptrOffset + posGlobalize(SA[leftTmp], limits)];

                const long tmp = rightTmp - leftTmp;
                if (tmp > 0)
                {
                    appendValue(nextBucketIndicesThread, Triple<TBktIndex, TBktIndex, unsigned>(leftTmp, rightTmp, newOffset));
                }
                leftTmp = rightTmp + 1;
            }

            leftTmp = middleValues.i2 + 1;
            while (leftTmp < end)
            {
                bptrOffset = getSeqNo(SA[leftTmp]) * bptrExtPerString;
                const TBktIndex rightTmp = bptr[bptrOffset + posGlobalize(SA[leftTmp], limits)];
                const long tmp = rightTmp - leftTmp;
                if (tmp > 0)
                {
                    appendValue(nextBucketIndicesThread, Triple<TBktIndex, TBktIndex, unsigned>(leftTmp, rightTmp, newOffset));
                }
                leftTmp = rightTmp + 1;
            }

            if (middleValues.i2 > middleValues.i1 + 1)
            {
                appendValue(nextBucketIndicesThread, Triple<TBktIndex, TBktIndex, unsigned>(middleValues.i1 + 1, middleValues.i2, 2*offset));
            }
        }

        //swap the oldBucketIndices with the new string for the next iteration
        clear(resultValues);
        clear(bucketIndicesSet);
        swap(bucketIndicesSet, nextBucketIndices);
        bucketIndices = concat(bucketIndicesSet);
    }

#if (SEQAN_ENABLE_DEBUG || SEQAN_ENABLE_TESTING) && SEQAN_ENABLE_PARALLELISM
        std::cout << " bucketSort: " << float(omp_get_wtime() - begin_time) << std::endl;

        const clock_t begin_time2 = omp_get_wtime();
#endif

    /*
     * We now have sorted all Buckets, except those starting with the same first and second character.
     * Theese Buckets can be sorted by using the order of the other buckets.
     * (This is a kind of seard copy)
     */

    //Now use seward copy to generate the order for all buckets starting with the same two characters
    SEQAN_OMP_PRAGMA(parallel for)
    for (unsigned i = 0; i < ALPHABETSIZE; ++i)
    {
        const TAlpha firstChar = 1 + i;

        TBktIndex leftVal = bkt[firstChar * bucketsInL1Bucket + firstChar * bucketsInL2Bucket];
        TBktIndex rightVal = bkt[firstChar * bucketsInL1Bucket + (firstChar + 1) * bucketsInL2Bucket];

        TBktIndex left = bkt[firstChar * bucketsInL1Bucket];
        TBktIndex right = bkt[(firstChar + 1) * bucketsInL1Bucket];

        //TODO seward kommentieren

        while (left < leftVal)
        {
            TSAValue tmp = SA[left]; //firstChar-bucket starts at this position of the text
            TTextSize seqOffset = getSeqOffset(tmp);
            const unsigned seqNum = getSeqNo(tmp);
            if (seqOffset <= 0)
            {
                ++left;
                continue;
            }
            const TAlpha character = 1 + ordValue(getSequenceByNo(seqNum, s)[--seqOffset]);
            if (firstChar == character)
            {
                const long bptrOffset = seqNum * bptrExtPerString;
                setSeqOffset(tmp, seqOffset);
                bptr[posGlobalize(tmp, limits) + bptrOffset] = leftVal;
                SA[leftVal] = tmp;
                ++leftVal;
            }
            ++left;
        }

        while (left < right)
        {
            --right;
            TSAValue tmp = SA[right];
            TTextSize seqOffset = getSeqOffset(tmp);
            const long seqNum = getSeqNo(tmp);
            TAlpha character;
            if (seqOffset > 0 && firstChar == (character = 1 + ordValue(getSequenceByNo(seqNum, s)[--seqOffset])))
            {
                --rightVal;
                const long bptrOffset = bptrExtPerString * seqNum;
                setSeqOffset(tmp, seqOffset);
                bptr[posGlobalize(tmp, limits) + bptrOffset] = rightVal;
                SA[rightVal] = tmp;
            }
        }
    }

#if (SEQAN_ENABLE_DEBUG || SEQAN_ENABLE_TESTING) && SEQAN_ENABLE_PARALLELISM
    std::cout << " seward: " << float(omp_get_wtime() - begin_time2)<<std::endl;
#endif
}

//update bptr after sort
//this part has its origin in the original algorithm by schürmann and stoye: updatePtrAndRefineBuckets_SaBucket
template<typename TSA, typename TBptr, typename TSize, typename TLimits>
Pair<TSize, TSize> refineBucket(TSA const & SA, TBptr const & bptr, String<Pair<TSize, TSize> > & nextBptr, TLimits const & limits,
        TSize const & left, TSize const & right, unsigned int const & offset, unsigned short const & d)
{
    const unsigned bptrExtPerString = (2 * d + 1);

    if (right - left == 1)
    {
        const unsigned bptrOffset1 = getSeqNo(SA[left]) * bptrExtPerString;
        appendValue(nextBptr, Pair<TSize, TSize>(bptrOffset1 + posGlobalize(SA[left], limits), left));
        return Pair<TSize, TSize>(0, 0);
    }

    typedef typename Value<TSA>::Type TSAVal;

    /*
     * from right to left compare sort keys, as long as they are equal they are in the same bucket.
     * set the bucket to the rightmost position
     */
    TSize leftInterval = right;
    TSize rightInterval = right;
    TSize tmp;

    while (left <= leftInterval
            && right < (tmp = getBptrVal(SA, bptr, limits, bptrExtPerString, offset, leftInterval)))// bptr[SA[leftInterval] + offset]))
    {
        do
        {
            unsigned bptrOffset = getSeqNo(SA[leftInterval]) * bptrExtPerString;
            // nextBptr[bptrOffset + posGlobalize(SA[leftInterval], limits)] =  rightInterval;
            appendValue(nextBptr, Pair<TSize, TSize>(bptrOffset + posGlobalize(SA[leftInterval], limits), rightInterval));
            --leftInterval;
        }
        while (left <= leftInterval && (TSize)getBptrVal(SA, bptr, limits, bptrExtPerString, offset, leftInterval)/*bptr[SA[leftInterval] + offset]*/
                == tmp);
        rightInterval = leftInterval;
    }

    rightInterval = leftInterval;
    while (left <= leftInterval && left <= (TSize)getBptrVal(SA, bptr, limits, bptrExtPerString, offset, leftInterval) //bptr[SA[leftInterval] + offset]
            && (TSize)getBptrVal(SA, bptr, limits, bptrExtPerString, offset, leftInterval) <= right)
    {
        unsigned bptrOffset = getSeqNo(SA[leftInterval]) * bptrExtPerString;
        //nextBptr[bptrOffset + posGlobalize(SA[leftInterval], limits)] =   rightInterval;
        appendValue(nextBptr, Pair<TSize, TSize>(bptrOffset + posGlobalize(SA[leftInterval], limits), rightInterval));
        --leftInterval;
    }

    const TSize middleRight = rightInterval;
    const TSize middleLeft = leftInterval;

    rightInterval = leftInterval;
    while (left <= leftInterval)
    {
        const TSize tmp2 = getBptrVal(SA, bptr, limits, bptrExtPerString, offset, leftInterval); //bptr[SA[leftInterval] + offset];
        do
        {
            unsigned bptrOffset = getSeqNo(SA[leftInterval]) * bptrExtPerString;
            //nextBptr[bptrOffset + posGlobalize(SA[leftInterval], limits)] =   rightInterval;
            appendValue(nextBptr,
                    Pair<TSize, TSize>(bptrOffset + posGlobalize(SA[leftInterval], limits), rightInterval));
            --leftInterval;
        }
        while (left <= leftInterval && (TSize)getBptrVal(SA, bptr, limits, bptrExtPerString, offset, leftInterval)/*bptr[SA[leftInterval] + offset] */
            == tmp2);
        rightInterval = leftInterval;
    }

    return Pair<TSize, TSize>(middleLeft, middleRight);
}

//finds the longest common prefix and increases the current offset
template<typename TSa, typename TBptr, typename TSize, typename TLimits>
unsigned int computeLCPAndOffset(TSa const & SA, TBptr const & bptr, TLimits const & limits, TSize const & left,
        TSize const & right, unsigned int const & offset, unsigned short const & d)
{
    typedef typename Value<TBptr>::Type TBptrVal;

    const unsigned bptrExtPerString = (2 * d + 1);

    unsigned int lcp = offset;
    while (true)
    {
        TSize index = left;
        unsigned bptrOffset = getSeqNo(SA[right]) * bptrExtPerString;
        const TBptrVal tmp = bptr[bptrOffset + posGlobalize(SA[right], limits) + lcp];
        while (index < right)
        {
            bptrOffset = getSeqNo(SA[index]) * bptrExtPerString;
            if (bptr[bptrOffset + posGlobalize(SA[index], limits) + lcp] != tmp)
            {
                return lcp;
            }
            ++index;
        }
        lcp += d;
    }
}

//Sort and refine a bucket with size 2 without recursion
template<typename TSA, typename TBptr, typename TIndex, typename TLimits>
void sortSizeTwo(TSA &SA, TBptr const & bptr, TLimits const & limits, unsigned short const & d, unsigned int const & offset,
        TIndex const & left, TIndex const & right)
{
    typedef typename Value<TSA>::Type TSAVal;

    const unsigned bptrExtPerString = (2 * d + 1);
    const unsigned bptrOffset1 = getSeqNo(SA[left]) * bptrExtPerString;
    const unsigned bptrOffset2 = getSeqNo(SA[right]) * bptrExtPerString;

    unsigned suffix1 = posGlobalize(SA[left], limits) + offset;
    unsigned suffix2 = posGlobalize(SA[right], limits) + offset;

    //find first non matching entry for both suffixes
    while (bptr[bptrOffset1 + suffix1] == bptr[bptrOffset2 + suffix2])
    {
        suffix1 += d;
        suffix2 += d;
    }

    //if necessary swap both values
    if (bptr[bptrOffset1 + suffix1] > bptr[bptrOffset2 + suffix2])
    {
        std::swap(SA[left], SA[right]);
    }
}

// gets the correct value from bptr, respects all necessary offsets
//
template<typename TSA, typename TBptrVal, typename TLimits, typename TExtString, typename TOffset,
        typename TIndex>
TBptrVal getBptrVal(TSA const & SA, String<TBptrVal> const & bptr, TLimits const & limits,
        TExtString const & bptrExtPerString, TOffset const & offset, TIndex const & index)
{
    //bptrOffset = getSeqNo(SA[index]) * bptrExtPerString
    return bptr[posGlobalize(SA[index], limits) + offset + getSeqNo(SA[index]) * bptrExtPerString];
}


// ============================================================================
// Functions
// ============================================================================

}  // namespace seqan

#endif  // #ifndef CORE_INCLUDE_SEQAN_INDEX_INDEX_SA_BPR_H_
