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

    //for each bucket, this contains it's startindex in SA
    String<TTextSize> bkt;

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

    sortAndRefinePhase2(SA, bptr, s, bkt, limits, d, getAlphabetSize(s)+1);

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

        //If the input has multiple sequences, there has to be an offset between each sequence
        //for a single input this is always 0
        long bptrOffset;

        //every thread has it's own bucketTable:
        TBkt &threadBucket = bktPerThread[threadNum];
        resize(threadBucket, bucketCount + 1, 0, Exact());

        //////////////////////////////
        // first iteration: compute rank (code_d) for each d-suffix between start and end and count the bucketsizes
        //////////////////////////////

        posLocalize(saValue, start, limits);
        bptrOffset = getSeqNo(saValue) * bptrExtPerString;

        //compute rank
        unsigned codeD = code_d(
                getSequenceByNo(getSeqNo(saValue), s), d, saValue,
                alphabeSizeWithDollar);
        bptr[start + bptrOffset] = codeD;//store rank
        threadBucket[codeD]++;//increase bucketSize

        for (TTextSize i = start + 1; i < end; ++i)
        {
            posLocalize(saValue, i, limits);
            bptrOffset = getSeqNo(saValue) * bptrExtPerString;

            //compute rank using previous value
            codeD = code_d(
                    getSequenceByNo(getSeqNo(saValue), s), d, saValue,
                    alphabeSizeWithDollar, tmpModulo,
                    codeD);
            bptr[i + bptrOffset] = codeD;//store rank

            threadBucket[codeD]++;//increase bucketSize
        }

        SEQAN_OMP_PRAGMA(barrier)

        //summarize buckets: each thread needs to know how many entries the previous threads are going to insert into each bucket
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
        //now bkt contains the starting positions for every bucket (=number of entrys in smaller buckets)
        //and bktPerThread stores how many entries the previous threads are going to insert into the bucket

        //////////////////////
        // insert the buckets into the suffixarry
        //////////////////////
        for (TTextSize i = start; i < end; ++i)
        {
            posLocalize(saValue, i, limits);
            bptrOffset = getSeqNo(saValue) * bptrExtPerString;
            //every thread fills his own range, computed with the bucketSize and the size of previous threads
            TBktValue index = (bkt[bptr[i + bptrOffset] + 1]
                    - threadBucket[bptr[i + bptrOffset]]) - 1;

            threadBucket[bptr[i + bptrOffset]]++;//increase number of inserted entrys
            SA[index] = saValue;//insert the value into the suffixarray
        }

        //up to now we misused bptr to save the rank-value of a suffix
        //now bptr should be a pointer from suffix to his bucket
        for (TTextSize i = start; i < end; ++i)
        {
            posLocalize(saValue, i, limits);
            bptrOffset = getSeqNo(saValue) * bptrExtPerString;
            //each bucketpointer points to the buckets last index
            //this is the next buckets startposition -1:
            bptr[i + bptrOffset] = bkt[bptr[i + bptrOffset] + 1] - 1;
        }

        SEQAN_OMP_PRAGMA(barrier)

        //the last d-suffixes of each sequence gets special handling corresponding to the sorting of the $-signs
        //this needs to be done to ensure stable sorting: each suffix containing one or more $-signs get a special bucket
        //this allows to have an order between the different $: $_1<$_2<...<$_x
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
                    //compute the rank
                    lastValue = code_d(getSequenceByNo(seq, s), d, i,
                            alphabeSizeWithDollar, tmpModulo, lastValue);

                    //this suffix gets an own bucket
                    bptr[bptrOffset + seqStartIndex + i] = bkt[lastValue];
                    bkt[lastValue] = bkt[lastValue] + 1;
                }

                //insert negative values after each sequence
                //this ensures unique sortkeys in refinement phase which results in a stable sort order between multiple sequences
                TTextSize j = seq * 2 * d + seq;
                j *= -1;
                for (TTextSize i = seqEndIndex; i <= seqEndIndex + 2 * d; ++i)
                {
                    bptr[seqStartIndex + i + bptrOffset] = --j;
                }
            }

            //some bkt-values got modified, this has to be undone
            //reinsert correct bkt values (they are later needed for seward copy)
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
//previousCodeD: rank of previous suffix
template<typename TText, typename TPos, typename TAlphabetSize>
unsigned code_d(TText const & s, unsigned short const & d, TPos const & pos, TAlphabetSize const & alphabetSize, unsigned const & modulo,
        unsigned const & previousCodeD)
{
    unsigned local = getSeqOffset(pos);
    if (local == 0) //don't use previous rank, if new sequence starts
        return code_d(s, d, pos, alphabetSize);

    if (length(s) <= local + d - 1)
        return alphabetSize * (previousCodeD % modulo) + 0;
    else
        return alphabetSize * (previousCodeD % modulo) + 1 + ordValue(s[local + d - 1]);
}

//Phase 2: sort and refine all Buckets
template<typename TSA, typename TText, typename TBptr, typename TBkt, typename TLimit, typename TD, typename TAlpha>
inline void sortAndRefinePhase2(TSA & SA, TBptr & bptr, TText const & s, TBkt const & bkt, TLimit const & limits, TD const d,
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

    ////////////////////////////
    //Initial: find all non empty Buckets and store their indices
    ////////////////////////////
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

    //////////////////////////////
    // Now sort and refine these Buckets
    //////////////////////////////

    StringSet<String<Triple<TBktIndex, TBktIndex, unsigned> > > nextBucketIndices;
    String<Triple<TBktIndex, TBktIndex, unsigned> > bucketIndices = concat(bucketIndicesSet);
    String<unsigned> split;
    String<Pair<long, long> > middleValues;//we need signed values

    String<String<Pair<TBktIndex, TBktIndex> > > nextBptrValuesPerThread;
    resize(nextBptrValuesPerThread, omp_get_max_threads());

    //minimal bucketcount of each job
    const long minBuckets = 100 * omp_get_max_threads();
    typedef _bprComparator<TTextSize, TSA, TBptr, TLimit, TTextSize> TSortFunctor;

    //sort as long as there are buckets to sort
    while (length(bucketIndices) > 0)
    {
        ///////////////////////
        // Split buckets in reasonable chunks,
        // each chunk gets sorted and then refined
        // (if the chunks are to big: refined values can't get used,
        //  if they are to small: to much parallel sorting can be used)
        ///////////////////////

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
        resize(middleValues, length(bucketIndices));

        ////////////////////////////
        //sort each split and refine afterwards in parallel
        //we can not sort and refine in parallel
        //if we wait to long with the refinement, the sorting gets ineffective.
        //so we use split to regularly switch between sorting and refinement.
        ////////////////////////////
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

                    //refine buckets
                    //middleValues contains the indices of those buckets, whose sortkeys are also within this bucket
                    middleValues[k] = refineBucket(SA, bptr, nextBptrValues, limits, start, end, offset, d);
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

            const Pair<long, long> middleValue = middleValues[k];

            String<Triple<TBktIndex, TBktIndex, unsigned> > &nextBucketIndicesThread = nextBucketIndices[omp_get_thread_num()];

            //increase offset
            unsigned int newOffset = offset + d;
            if ((TBktIndex)bptr[bptrOffset + posGlobalize(SA[start], limits)] == end)
            {
                //further increase offset
                newOffset = computeLCPAndOffset(SA, bptr, limits, start, end, newOffset, d);
            }

            //search for subbucket-indices within 'start' and 'end' and add them to the list
            //ignore range between the middlevalues, it gets treated separately
            TBktIndex leftTmp = start;
            while ((signed)leftTmp < middleValue.i1)
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

            leftTmp = middleValue.i2 + 1;
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

            //the range between these values is the bucket whose sortkeys were within the original bucket
            //hence we now that not only the suffixes but also the sortkeys share a common prefix of length 'offset'
            //so we can double the offset
            if (middleValue.i2 > middleValue.i1 + 1)
            {
                appendValue(nextBucketIndicesThread, Triple<TBktIndex, TBktIndex, unsigned>(middleValue.i1 + 1, middleValue.i2, 2 * offset));
            }
        }

        //swap the oldBucketIndices with the new string for the next iteration
        clear(middleValues);
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
     * (This is based on sewards copy algorithm)
     */

    SEQAN_OMP_PRAGMA(parallel for)
    for (unsigned i = 0; i < ALPHABETSIZE; ++i)
    {
        const TAlpha firstChar = 1 + i;

        //this finds all buckets starting with the current character twice (interval: [AA;AB[ )
        //(this are the buckets we want to sort)
        TBktIndex leftVal = bkt[firstChar * bucketsInL1Bucket + firstChar * bucketsInL2Bucket];
        TBktIndex rightVal = bkt[firstChar * bucketsInL1Bucket + (firstChar + 1) * bucketsInL2Bucket];

        //and here we search for all buckets starting with the current character (interval [A;B[ )
        //(this are all buckets we need to retrive the order of above buckets)
        TBktIndex left = bkt[firstChar * bucketsInL1Bucket];
        TBktIndex right = bkt[(firstChar + 1) * bucketsInL1Bucket];


        //search for the smallest (sortorder) character whose predecessor is the same as the character.
        //this has to be the smallest suffix with the character twice, and so on
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

            //check if previous character is the same
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

        //now do the same from the right side
        while (left < right)
        {
            --right;
            TSAValue tmp = SA[right];
            TTextSize seqOffset = getSeqOffset(tmp);
            const long seqNum = getSeqNo(tmp);
            TAlpha character;
            //check if previous character is the same
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
//this function is a modified version of updatePtrAndRefineBucketsNum_SaBucket written by
//schürmann and stoye for the original bpr implementation bpr-0.9.0
//returning pair describes the range of the subbucket whose sortkeys are also within this bucket
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
    long leftInterval = right;
    long rightInterval = right;
    TSize tmp;

    //at first only refine buckets whose sortkey (bptr+offset) is larger than 'right'
    while ((signed)left <= leftInterval
            && right < (tmp = getBptrVal(SA, bptr, limits, bptrExtPerString, offset, leftInterval)))// bptr[SA[leftInterval] + offset]))
    {
        do
        {
            unsigned bptrOffset = getSeqNo(SA[leftInterval]) * bptrExtPerString;
            // nextBptr[bptrOffset + posGlobalize(SA[leftInterval], limits)] =  rightInterval;
            appendValue(nextBptr, Pair<TSize, TSize>(bptrOffset + posGlobalize(SA[leftInterval], limits), rightInterval));
            --leftInterval;
        }
        while ((signed)left <= leftInterval && (TSize)getBptrVal(SA, bptr, limits, bptrExtPerString, offset, leftInterval)/*bptr[SA[leftInterval] + offset]*/
                == tmp);
        rightInterval = leftInterval;
    }

    //since bucketpointer between left and right could change in the other steps in this function (though the offset), they have to be treated separately
    rightInterval = leftInterval;
    while ((signed)left <= leftInterval && left <= (TSize)getBptrVal(SA, bptr, limits, bptrExtPerString, offset, leftInterval) //bptr[SA[leftInterval] + offset]
            && (TSize)getBptrVal(SA, bptr, limits, bptrExtPerString, offset, leftInterval) <= right)
    {
        unsigned bptrOffset = getSeqNo(SA[leftInterval]) * bptrExtPerString;
        //nextBptr[bptrOffset + posGlobalize(SA[leftInterval], limits)] =   rightInterval;
        appendValue(nextBptr, Pair<TSize, TSize>(bptrOffset + posGlobalize(SA[leftInterval], limits), rightInterval));
        --leftInterval;
    }

    const TSize middleRight = rightInterval;
    const TSize middleLeft = leftInterval;

    //now refine the remainig range
    rightInterval = leftInterval;
    while ((signed)left <= leftInterval)
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
        while ((signed)left <= leftInterval && (TSize)getBptrVal(SA, bptr, limits, bptrExtPerString, offset, leftInterval)/*bptr[SA[leftInterval] + offset] */
            == tmp2);
        rightInterval = leftInterval;
    }

    return Pair<long, long>(middleLeft, middleRight);
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
