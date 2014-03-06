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

#ifndef CORE_INCLUDE_SEQAN_INDEX_INDEX_BWT_BCR_H_
#define CORE_INCLUDE_SEQAN_INDEX_INDEX_BWT_BCR_H_

namespace seqan
{

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

template<typename TPos, typename InType, typename OutType>
struct _posComparator : public std::unary_function<InType, OutType>
{

    TPos const &s;

    _posComparator(TPos const &_s) :
            s(_s)
    {
    }

    OutType operator()(InType index)
    {
        return s[index].i2;
    }
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================


/**
.Function.createBwt:
..summary:Computes the BWT for given input
..description:
Creates the Burrows–Wheeler transform for the input text. The algorithm is based on the lightweight Algorithm BCR by Bauer, Cox and Rosone
( http://dx.doi.org/10.1016/j.tcs.2012.02.002 ).
It was modified to work completely in memory and make use of multiple cores.
..createBwt(BWT, sentinelPosition, text, sentinelSub)
..param.BWT: result BWT, it has to have the correct size an type
..param.sentinelPosition: string that stores the positions of the sentinels in the bwt.
..param.text: the Inputtext
..param.sentinelSub: character which replaces the sentinel char in the bwt
*/
template<typename TBWT, typename TText, typename TSentinelPosition, typename TSentinelSub>
void createBwt(TBWT & BWT, TSentinelPosition & sentinelPos, TText & text, TSentinelSub const sentinelSub)
{
    ModifiedString<TText> shallowCopy(text);
    StringSet<ModifiedString<TText> > set;
    appendValue(set, shallowCopy);

    createBwt(BWT, sentinelPos, set, sentinelSub);
}

/**
.Function.createBwt:
..summary:Computes the BWT for given input
..description:
Creates the Burrows–Wheeler transform for the input text. The algorithm is based on the lightweight Algorithm BCR by Bauer, Cox and Rosone
( http://dx.doi.org/10.1016/j.tcs.2012.02.002 ).
It was modified to work completely in memory and make use of multiple cores.
..createBwt(BWT, sentinelPosition, text, sentinelSub)
..param.BWT: result BWT, it has to have the correct size an type
..param.sentinelPosition: string that stores the positions of the sentinels in the bwt.
..param.text: the Inputtexts
..param.sentinelSub: character which replaces the sentinel char in the bwt
*/
template<typename TBWT, typename TText, typename TSentinelPosition, typename TSentinelSub>
void createBwt(TBWT & BWT, TSentinelPosition & sentinelPos, StringSet<TText> & text, TSentinelSub const sentinelSub)
{
    typedef typename Value<TBWT>::Type TAlphabet;
    typedef typename ValueSize<TAlphabet>::Type TAlphabetSize;
    typedef typename Size<TText>::Type TValueSize;

    TAlphabetSize const ALPHABETSIZE = ValueSize<TAlphabet>::VALUE;
    TAlphabetSize const ALPHABETSIZEWITHDOLLAR = ALPHABETSIZE + 1;
    TValueSize const count = countSequences(text);

    resize(sentinelPos, count, Exact());

    //get longest sequence
    TValueSize maxLength = 0;
    for (unsigned i = 0; i < length(text); ++i)
    {
        maxLength = _max(maxLength, length(text[i]));
    }

    unsigned bs = 1;
    unsigned bc = ALPHABETSIZEWITHDOLLAR;

#if SEQAN_ENABLE_PARALLELISM
    //Default: 1 bucket per character
    //With more Threads: create buckets for pairs, triples or more
    unsigned threadCount = omp_get_max_threads();
    unsigned minBuckets = 200 * threadCount;
    while (bc < minBuckets)
    {
        bs += 1;
        bc *= ALPHABETSIZEWITHDOLLAR;
    }
#endif

    unsigned const bucketSize = bs;
    unsigned const bucketCount = bc;

    //Triple:
    // 1: index for sorting the entries (should match index in String)
    // 2: the character in the BWT
    // 3: true, if the character is a sentinel
    typedef Triple<unsigned, TAlphabet, bool> BucketValueType;

    //TODO: replace short/int/long by Size<>::Type / Value<>...

    String<String<BucketValueType> > bucket;
    resize(bucket, bucketCount, Exact());
    String<unsigned> bucketVolume; //current size of each bucket
    resize(bucketVolume, bucketCount, 0);

    // for each bucket: string with length count * alphabetsize which contains for each position the count of each character up to this position
    String<String<unsigned> > charCountBuffer;
    resize(charCountBuffer, bucketCount);
    String<unsigned> charCountBufferSize; //current buffer size
    resize(charCountBufferSize, bucketCount, 0);

    //Allocate the bucketmemory here, to avoid append and resize operations in parallel sections
    preallocateMemory(text, bucketCount, bucketSize, ALPHABETSIZE, bucket, charCountBuffer);

    /*
     * P: Position in Bucket to Insert
     * N: Index of String in (input)text
     * Q: Number of Bucket to Insert
     */
    String<Pair<TValueSize, TValueSize> > pos; //N, P
    resize(pos, count, Exact());
    String<TValueSize> qIndex;
    resize(qIndex, count, Exact()); //Q

    String<TValueSize> tempQIndex; //Q in next Iteration
    resize(tempQIndex, count, Exact());

    //number of inserts in last iteration for each Bucket
    //since we allocated all memory, we can not use the string-length
    String<unsigned> bucketInsertCount;
    resize(bucketInsertCount, bucketCount, 0, Exact());

    /*
     For each iteration build:
     B_j(h): bwt for characters with suffix starting with 'h'
     N_j(h): Array containing Index i of S of the j-1'th suffix
     P_j(h): Array contains the position of j-1'th suffix in B_j(h)
     */

    TValueSize iteration = 0;

    TValueSize i;
    SEQAN_OMP_PRAGMA(parallel for)
    for (i = 0; i < count; ++i)
    {
        /*
         * This defines the weight of the $-signs
         * tIndex = i: Sort $_1 < $_2 < ... < $_n
         * tIndex = count - 1 - i: Sort $_1 > $_2 > ... > $_n
         */
        const TValueSize tIndex = count - 1 - i;

        const TText& currentText = text[tIndex];
        const bool isLastCharOfText = iteration == length(currentText);
        const TAlphabet newChar = (isLastCharOfText) ? sentinelSub : currentText[length(currentText) - iteration - 1];

        //insert last char of text
        bucket[0][i] = BucketValueType(i, newChar, isLastCharOfText);

        pos[i].i1 = tIndex;
        pos[i].i2 = i;
        if (isLastCharOfText)
        {
            //set invalid bucketnummer to ignore in further iterations
            qIndex[i] = bucketCount;
            tempQIndex[i] = bucketCount;
        }
        else
        {
            qIndex[i] = 0;
        }
    }
    bucketVolume[0] = count;
    bucketInsertCount[0] = count;

#if (SEQAN_ENABLE_DEBUG || SEQAN_ENABLE_TESTING) && SEQAN_ENABLE_PARALLELISM
    double countTime = 0;
    double insertTime = 0;
    double getPosTime = 0;
    double getSortTime = 0;
    double getSort2Time = 0;
#endif

    //Store Insertcoint per thread to avoid locking
    String<String<unsigned> > tempBucketInsertCount;
    resize(tempBucketInsertCount, omp_get_max_threads(), Exact());
    for (unsigned int i = 0; i < length(tempBucketInsertCount); ++i)
    {
        resize(tempBucketInsertCount[i], bucketCount, Exact());
    }

    //Buffer values while sorting newly inserted
    String<String<BucketValueType> > buffer;
    resize(buffer, bucketCount, Exact());

    typedef _posComparator<String<Pair<TValueSize, TValueSize> >, unsigned, TValueSize> TPosSortFunctor;
    TPosSortFunctor sortFunctor = TPosSortFunctor(pos);

    const unsigned powAlphabetBucketSize = pow(ALPHABETSIZEWITHDOLLAR, bucketSize - 1);

    //Itearate characterwise and insert them into the buckets
    for (iteration = 1; iteration <= maxLength; ++iteration)
    {

    SEQAN_OMP_PRAGMA(parallel shared(text, bucket, bucketVolume, charCountBuffer, charCountBufferSize, tempBucketInsertCount, bucketInsertCount, pos, qIndex, tempQIndex, buffer))
    {
        String<unsigned> &threadBktInsertCount = tempBucketInsertCount[omp_get_thread_num()];

        //reset bucketInsertCounter
        for (long bucketIndex = 0; bucketIndex < bucketCount; ++bucketIndex)
        {
            threadBktInsertCount[bucketIndex] = 0;
        }

#if (SEQAN_ENABLE_DEBUG || SEQAN_ENABLE_TESTING) && SEQAN_ENABLE_PARALLELISM
        double begin_count_time = 0;
        if(omp_get_thread_num()==0)
        {
            begin_count_time = omp_get_wtime();
        }
#endif

        //count characters per bucket for each position and store count-values in charCountBuffer
        SEQAN_OMP_PRAGMA(for schedule(dynamic))
        for (long bucketIndex = 0; bucketIndex < bucketCount; ++bucketIndex)
        {
            String<unsigned> &currentBuffer = charCountBuffer[bucketIndex];
            String<BucketValueType> &currentBucket = bucket[bucketIndex];

            const int oldbucketLength = charCountBufferSize[bucketIndex];
            const int bucketLength = bucketVolume[bucketIndex];

            if(oldbucketLength == bucketLength)
            {
                //this bucket did not change, no need to recount
                continue;
            }
            charCountBufferSize[bucketIndex] = bucketLength;

            //Get the first changed index since last iteration. We don't need to recount the previous values
            const unsigned startBucketIndex = bucketIndex == 0 ? 0 : bucketInsertCount[bucketIndex - 1];
            const TValueSize &currentPos = pos[startBucketIndex].i2;

            //Initially all counts are 0.
            if (currentPos == 0)
            {
                for (unsigned j = 0; j < ALPHABETSIZE; ++j)
                {
                    currentBuffer[j] = 0;
                }
                BucketValueType &bktValue = currentBucket[0];
                if (!bktValue.i3)
                {
                    ++currentBuffer[ordValue(bktValue.i2)];
                }
            }


            for (int index = (currentPos==0 ? 1 : currentPos); index < bucketLength; ++index)
            {
                unsigned currentBufferStart = index * ALPHABETSIZE;
                unsigned lastBufferStart = currentBufferStart - ALPHABETSIZE;

                //First copy the values from previous entry
                for (unsigned j = 0; j < ALPHABETSIZE; ++j)
                {
                    currentBuffer[currentBufferStart + j] = currentBuffer[lastBufferStart + j];
                }

                //Increment count for current character
                BucketValueType &bktValue = currentBucket[index];
                if (!bktValue.i3)
                {
                    ++currentBuffer[currentBufferStart + ordValue(bktValue.i2)];
                }
            }
        }

#if (SEQAN_ENABLE_DEBUG || SEQAN_ENABLE_TESTING) && SEQAN_ENABLE_PARALLELISM
        if(omp_get_thread_num()==0)
        {
            countTime += omp_get_wtime() - begin_count_time;
        }
        double begin_getpos_time = 0;
        if(omp_get_thread_num()==0)
        {
            begin_getpos_time = omp_get_wtime();
        }
#endif

        //get insert positions for the next character in every input-text
        //Do this in parallel for every bucket
        SEQAN_OMP_PRAGMA(for schedule(guided))
        for (unsigned bucketIndex = 0; bucketIndex < bucketCount; ++bucketIndex)
        {
            const unsigned startBucketIndex = bucketIndex == 0 ? 0 : bucketInsertCount[bucketIndex - 1];
            const unsigned endBucketIndex = bucketInsertCount[bucketIndex];

            const String<BucketValueType> &currentBucket = bucket[bucketIndex];

            //all positions corresponding to this Bucket
            for (unsigned i = startBucketIndex; i < endBucketIndex; ++i)
            {
                if (bucketIndex != qIndex[i])
                {
                    //this happens only if the text ended in the previous iteration
                    //It will be ignored in further iterations
                    continue;
                }

                const TValueSize &currentPos = pos[i].i2;
                const BucketValueType &bv = currentBucket[currentPos]; //in previous iteration inserted char
                const unsigned charValue = ordValue(bv.i2);

                //count occurences in every bucket until we reach current position
                //(For bucketsize>0 it is sufficient to count only the buckets starting with the same character)
                const long startBucket = bucketIndex - (bucketIndex % ALPHABETSIZEWITHDOLLAR);
                long charCount = 0;
                for (unsigned j = startBucket; j < bucketIndex; ++j)
                {
                    String<unsigned> &cur = charCountBuffer[j];
                    const unsigned curBufferSize = charCountBufferSize[j];
                    if (curBufferSize > 0)
                    {
                        charCount += cur[ALPHABETSIZE * (curBufferSize - 1) + charValue];
                    }
                }
                if (currentPos > 0)
                {
                    charCount += charCountBuffer[bucketIndex][ALPHABETSIZE * (currentPos - 1) + charValue];
                }

                //remove last char of bucket (shift position (division)) and append new char in front (pow)
                const long newBucketIndex = (startBucket / ALPHABETSIZEWITHDOLLAR)
                        + (powAlphabetBucketSize * (1 + charValue));//1+ to ignore the firstBucket

                pos[i].i2 = charCount;
                tempQIndex[i] = newBucketIndex;
                ++threadBktInsertCount[newBucketIndex];
            }
        }

        //accumulate bucket insert count of every thread
        SEQAN_OMP_PRAGMA(for)
        for (long bucketIndex = 0; bucketIndex < bucketCount; ++bucketIndex)
        {
            for (int t = 1; t < omp_get_num_threads(); ++t)
            {
                tempBucketInsertCount[0][bucketIndex] += tempBucketInsertCount[t][bucketIndex];
            }
        }

        //set bucket insert count
        SEQAN_OMP_PRAGMA(single)
        {
            bucketInsertCount[0] = tempBucketInsertCount[0][0];
            for (long bucketIndex = 1; bucketIndex < bucketCount; ++bucketIndex)
            {
                bucketInsertCount[bucketIndex] = bucketInsertCount[bucketIndex - 1]
                + tempBucketInsertCount[0][bucketIndex];
            }
        }

        //resize Buffer, to avoid resize in parallel section
        SEQAN_OMP_PRAGMA(single nowait )
        {
            for (long bucketIndex = 0; bucketIndex < bucketCount; ++bucketIndex)
            {
                const unsigned startBucketIndex = bucketIndex == 0 ?
                        0 : bucketInsertCount[bucketIndex - 1];
                const unsigned endBucketIndex = bucketInsertCount[bucketIndex];
                resize(buffer[bucketIndex], endBucketIndex - startBucketIndex);
            }
        }

#if (SEQAN_ENABLE_DEBUG || SEQAN_ENABLE_TESTING) && SEQAN_ENABLE_PARALLELISM
        if(omp_get_thread_num()==0)
        {
            getPosTime += omp_get_wtime() - begin_getpos_time;
        }
        double begin_sort_time = 0;
        if(omp_get_thread_num()==0)
        {
            begin_sort_time = omp_get_wtime();
        }
#endif

        //Swap current and next qIndex, then sort qIndex, tempQIndex and pos according to the values in qIndex
        SEQAN_OMP_PRAGMA(single nowait)
        {
            swap(qIndex, tempQIndex);
            quickSort(qIndex, tempQIndex, pos, 0UL, length(qIndex) - 1, 0);
        }

        //Wait for parallel Quicksort tasks to finish
        SEQAN_OMP_PRAGMA(barrier)

#if (SEQAN_ENABLE_DEBUG || SEQAN_ENABLE_TESTING) && SEQAN_ENABLE_PARALLELISM
        if(omp_get_thread_num()==0)
        {
            getSortTime += omp_get_wtime() - begin_sort_time;
        }
        double begin_sort2_time = 0;
        if(omp_get_thread_num()==0)
        {
            begin_sort2_time = omp_get_wtime();
        }
#endif

        //We need to insert the new values per Bucket in ascending order.
        SEQAN_OMP_PRAGMA(for nowait)
        for (long bucketIndex = 0; bucketIndex < bucketCount; ++bucketIndex)
        {
            //get new position in ascending order
            const unsigned startBucketIndex = bucketIndex == 0 ?
                    0 : bucketInsertCount[bucketIndex - 1];
            const unsigned endBucketIndex = bucketInsertCount[bucketIndex];
            if (endBucketIndex > startBucketIndex)
                //sorts the buckets section in pos
                doQuickSort(qsortParallel(), pos, sortFunctor, startBucketIndex, endBucketIndex - 1);
        }

        //Wait for Quicksort tasks to finish
        SEQAN_OMP_PRAGMA(barrier)

#if (SEQAN_ENABLE_DEBUG || SEQAN_ENABLE_TESTING) && SEQAN_ENABLE_PARALLELISM
        if(omp_get_thread_num()==0)
        {
            getSort2Time += omp_get_wtime() - begin_sort2_time;
        }

        double begin_insert_time = 0;
        if(omp_get_thread_num()==0)
        {
            begin_insert_time = omp_get_wtime();
        }
#endif

        //Now insert the new characters into the buckets, in parallel for each bucket
        SEQAN_OMP_PRAGMA(for schedule(guided))
        for (long bucketIndex = 0; bucketIndex < bucketCount; ++bucketIndex)
        {
            const unsigned startBucketIndex = bucketIndex == 0 ?
                    0 : bucketInsertCount[bucketIndex - 1];
            const unsigned endBucketIndex = bucketInsertCount[bucketIndex];

            String<BucketValueType> &currentBuffer = buffer[bucketIndex];
            String<BucketValueType> &currentBucket = bucket[bucketIndex];
            unsigned &currentBucketVolume = bucketVolume[bucketIndex];

            for (unsigned i = startBucketIndex; i < endBucketIndex; ++i)
            {
                const TText& currentText = text[pos[i].i1];
                const bool isLastCharOfText = iteration == length(currentText);

                //replace with sentinelSub
                const TAlphabet newChar = (isLastCharOfText) ? sentinelSub :
                        currentText[length(currentText) - (int) iteration - 1];

                //Insert new values into the buffer. This avoids many move operations within the array
                BucketValueType &cur = currentBuffer[i - startBucketIndex];

                cur.i1 = pos[i].i2;
                cur.i2 = newChar;
                cur.i3 = isLastCharOfText;

                if (isLastCharOfText)
                {
                    //invalidate to ignore this input-text index in next iterations
                    qIndex[i] = bucketCount;
                    tempQIndex[i] = bucketCount;
                }
            }

            //now sort the new values from buffer to the correct positions
            sortBwtBucket(currentBucket, currentBucketVolume, currentBuffer, endBucketIndex - startBucketIndex);
            currentBucketVolume += endBucketIndex - startBucketIndex;
        }

#if (SEQAN_ENABLE_DEBUG || SEQAN_ENABLE_TESTING) && SEQAN_ENABLE_PARALLELISM
        if(omp_get_thread_num()==0)
        {
            insertTime += omp_get_wtime() - begin_insert_time;
        }
#endif

    }
    }

#if (SEQAN_ENABLE_DEBUG || SEQAN_ENABLE_TESTING) && SEQAN_ENABLE_PARALLELISM
    std::cout << " countTime:   " << countTime << std::endl;
    std::cout << " getPosTime:   " << getPosTime << std::endl;
    std::cout << " getSortTime:   " << getSortTime << std::endl;
    std::cout << " getSort2Time:   " << getSort2Time << std::endl;
    std::cout << " insertValue:   " << insertTime << std::endl;
#endif

    //write all values to the output-string
    int index = 0;
    int sentinelIndex = 0;
    for (long bucketIndex = 0; bucketIndex < bucketCount; ++bucketIndex)
    {
        String<BucketValueType> &currentBucket = bucket[bucketIndex];
        for (unsigned j = 0; j < length(currentBucket); ++j)
        {

            BucketValueType &bv = currentBucket[j];
            if (bv.i3)
            {
                sentinelPos[sentinelIndex++] = index;
            }
            BWT[index++] = bv.i2;
        }
    }
}

//allocate maximum size for each needed string. This avoids many resizes especially in parallel sections
template<typename TText, typename TBucketCount, typename TBucketSize, typename TAlphabetSize, typename TBucket,
    typename TCharCountBuffer>
void preallocateMemory(TText const & text, TBucketCount const & bucketCount, TBucketSize const & bucketSize,
        TAlphabetSize const & alphabetSize, TBucket & bucket, TCharCountBuffer & charCountBuffer)
{
    typedef typename Size<TText>::Type TValueSize;
    TAlphabetSize const alphabetSizeWithDollar = alphabetSize + 1;

#if (SEQAN_ENABLE_DEBUG || SEQAN_ENABLE_TESTING) && SEQAN_ENABLE_PARALLELISM
    double beginAllocTime = omp_get_wtime();
#endif

    String<unsigned> bucketCapacity;
    resize(bucketCapacity, bucketCount, 0);

    //Compute Bucketsize in parrallel
    TValueSize sequence;
    SEQAN_OMP_PRAGMA(parallel for)
    for(sequence = 0; sequence < countSequences(text); ++sequence)
    {
        for(unsigned t = 0; t < length(text[sequence]); ++t)
        {
            unsigned charValue = 0;
            for(TBucketSize b = 0; b < bucketSize; ++b)
            {
                if(length(text[sequence]) > t + b)
                charValue += (1 + ordValue(text[sequence][t + b])) * pow(alphabetSizeWithDollar, bucketSize - b - 1);
            }

            SEQAN_OMP_PRAGMA(atomic update)
            ++bucketCapacity[charValue];
        }
    }
    bucketCapacity[0] = countSequences(text);

    //allocate memory sequential, since it is very slow with openMp
    for (unsigned s = 0; s < bucketCount; ++s)
    {
        resize(bucket[s], bucketCapacity[s], Exact());
        resize(charCountBuffer[s], bucketCapacity[s] * alphabetSize, Exact());
    }

#if (SEQAN_ENABLE_DEBUG || SEQAN_ENABLE_TESTING) && SEQAN_ENABLE_PARALLELISM
    std::cout << " AllocTime:   " << omp_get_wtime() - beginAllocTime << std::endl;
#endif
}

/*
 * We have 'insertCount' new entrys in bucket. Sort them after with their position in .i1
 */
template<typename BucketValueType>
void sortBwtBucket(String<BucketValueType> & bucket, unsigned const & bucketVolume,
        String<BucketValueType> const & buffer, unsigned const & bufferVolume)
{
    if (bufferVolume == 0)
        return;

    unsigned insertCount = 0;

    int curBucketIndex = bucketVolume - 1;
    int curBufferIndex = bufferVolume - 1;
    int curWriteIndex = bucketVolume + bufferVolume - 1;

    //do a kind of merge sort to insert buffer into bucket, starting with the last values
    while (bufferVolume > insertCount)
    {
        BucketValueType const & curBuf = buffer[curBufferIndex];

        if (curBucketIndex >= 0)
        {
            BucketValueType &curBuk = bucket[curBucketIndex];

            if ((int) (curBuf.i1 - curBufferIndex) <= (int) curBuk.i1)
            {
                //write bucket value
                bucket[curWriteIndex] = curBuk;
                --curBucketIndex;

            }
            else
            {
                //write buffer value
                bucket[curWriteIndex] = curBuf;
                --curBufferIndex;
                ++insertCount;
            }
        }
        else
        {
            bucket[curWriteIndex] = curBuf;
            --curBufferIndex;
            ++insertCount;
        }
        bucket[curWriteIndex].i1 = curWriteIndex;
        --curWriteIndex;
    }
}

//do insertion sort for qIndex, move values from qIndexTemp and pos accordingly
template<typename TAlphabet, typename TFoo1, typename TFoo2, typename TIndex>
inline void insertionSort(String<TAlphabet> & qIndex, String<TFoo1> & qIndexTemp, String<TFoo2> & pos,
        TIndex const & left, TIndex const & right)
{
    TAlphabet keyQ;
    TFoo1 keyQTemp;
    TFoo2 keyPos;
    TIndex j;

    for (TIndex i = left + 1; i <= right; ++i)
    {
        keyQ = qIndex[i];
        keyQTemp = qIndexTemp[i];
        keyPos = pos[i];

        j = i - 1;

        int c = 1;
        while ((j >= left) && (qIndex[j] > keyQ))
        {
            qIndex[j + 1] = qIndex[j];
            qIndexTemp[j + 1] = qIndexTemp[j];
            pos[j + 1] = pos[j];
            if (j == 0)
            {
                c = 0;
                break; //underflow with unsigned
            }
            --j;
        }

        qIndex[j + c] = keyQ;
        qIndexTemp[j + c] = keyQTemp;
        pos[j + c] = keyPos;
    }
}

//do quicksort sort for qIndex, move values from qIndexTemp and pos accordingly
//create tasks until recursion is too deep. If length gets shorter than 15 use insertionSort instead
template<typename TAlphabet, typename TFoo1, typename TFoo2, typename TIndex>
void quickSort(String<TAlphabet> & qIndex, String<TFoo1> & qIndexTemp, String<TFoo2> & pos, TIndex const left,
        TIndex const right, long const depth)
{
    if (right - left < 15)
    {
        insertionSort(qIndex, qIndexTemp, pos, left, right);
        return;
    }

    const TIndex halfBucketSize = (right - left) / 2;
    TAlphabet pivot = qIndex[right];
    const TAlphabet pivotB = qIndex[left];
    const TAlphabet pivotC = qIndex[left + halfBucketSize];

    const int medianNumber = medianOfThree(pivot, pivotB, pivotC);

    if (medianNumber != 0)
    {
        pivot = medianNumber == 1 ? pivotB : pivotC;

        const TIndex swapIndex = medianNumber == 1 ? left : left + halfBucketSize;
        std::swap(qIndex[swapIndex], qIndex[right]);
        std::swap(qIndexTemp[swapIndex], qIndexTemp[right]);
        std::swap(pos[swapIndex], pos[right]);
    }

    //lomutos partitioning scheme
    TIndex leftTmp = left;
    TIndex rightTmp = left;
    while (rightTmp < right)
    {
        if (qIndex[rightTmp] <= pivot)
        {
            std::swap(qIndex[leftTmp], qIndex[rightTmp]);
            std::swap(qIndexTemp[leftTmp], qIndexTemp[rightTmp]);
            std::swap(pos[leftTmp], pos[rightTmp]);

            leftTmp++;
        }
        rightTmp++;
    }

    std::swap(qIndex[leftTmp], qIndex[rightTmp]);
    std::swap(qIndexTemp[leftTmp], qIndexTemp[rightTmp]);
    std::swap(pos[leftTmp], pos[rightTmp]);

    if (depth > 10)
    {
        //do sequential recursion
        if (right > leftTmp + 1)
        {
            quickSort(qIndex, qIndexTemp, pos, leftTmp + 1, right, depth);
        }
        do
        {
            --leftTmp;
        }
        while (leftTmp > left && qIndex[leftTmp] == pivot);
        if (leftTmp > left)
        {
            quickSort(qIndex, qIndexTemp, pos, left, leftTmp, depth);
        }
        return;
    }

    //create task for parallel recursion
    SEQAN_OMP_PRAGMA(task firstprivate(right, leftTmp, depth) shared(pos, qIndex, qIndexTemp))
    {
        if (right > leftTmp + 1)
        {
            quickSort(qIndex, qIndexTemp, pos, leftTmp + 1, right, depth +1);
        }
    }
    do
    {
        --leftTmp;
    }
    while (leftTmp > left && qIndex[leftTmp] == pivot);
    if (leftTmp > left)
    {
        quickSort(qIndex, qIndexTemp, pos, left, leftTmp, depth + 1);
    }
}

}  // namespace seqan

#endif  // #ifndef CORE_INCLUDE_SEQAN_INDEX_INDEX_BWT_BCR_H_
