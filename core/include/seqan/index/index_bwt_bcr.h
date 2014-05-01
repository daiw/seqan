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
    typedef Triple<TValueSize, TAlphabet, bool> BucketValueType;

    String<String<BucketValueType> > bucket;
    resize(bucket, bucketCount, Exact());
    String<TValueSize> bucketVolume; //current size of each bucket (since we allocate all memory at once, we can not use the string-length for this)
    resize(bucketVolume, bucketCount, 0);

    // for each bucket: string with length count * alphabetsize which contains for each position the count of each character up to this position
    String<String<TValueSize> > charCountBuffer;
    resize(charCountBuffer, bucketCount);
    String<TValueSize> charCountBufferSize; //current buffer size
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
    String<TValueSize> bucketInsertCount;
    resize(bucketInsertCount, bucketCount, 0, Exact());

    /*
     For each iteration build:
     B_j(h): bwt for characters with suffix starting with 'h' -> is in 'bucket'
     N_j(h): Array containing Index i of S of the j-1'th suffix -> is the first value in 'pos'
     P_j(h): Array contains the position of j-1'th suffix in B_j(h) -> is the second value in 'pos'
     */

    /*
     * The first iteration is easy: insert the last char of the text in bucket 0
     * (since they are followed by the sentinel character)
     */
    TValueSize iteration = 0;
    SEQAN_OMP_PRAGMA(parallel for)
    for (TValueSize i = 0; i < count; ++i)
    {
        /*
         * This defines the weight of the $-signs
         * tIndex = i: Sort $_1 < $_2 < ... < $_n
         * tIndex = count - 1 - i: Sort $_1 > $_2 > ... > $_n
         */
        TValueSize const tIndex = count - 1 - i;

        TText const & currentText = text[tIndex];
        bool const isLastCharOfText = iteration == length(currentText);
        TAlphabet const newChar = (isLastCharOfText) ? sentinelSub : currentText[length(currentText) - iteration - 1];

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

    //Store Insertcount per thread to avoid locking
    String<String<TValueSize> > tempBucketInsertCount;
    resize(tempBucketInsertCount, omp_get_max_threads(), Exact());
    typedef typename Size<String<String<TValueSize> > >::Type TBktInsertSize;
    for (TBktInsertSize i = 0; i < length(tempBucketInsertCount); ++i)
    {
        resize(tempBucketInsertCount[i], bucketCount, Exact());
    }

    String<String<BucketValueType> > buffer;//Buffer is needed for efficiently inserting multiple new values into the buckets
    resize(buffer, bucketCount, Exact());

    typedef _posComparator<String<Pair<TValueSize, TValueSize> >, TValueSize, TValueSize> TPosSortFunctor;
    TPosSortFunctor sortFunctor = TPosSortFunctor(pos); //sortfunctor, which allows to sort pos, by the insert position within the bucket

    const unsigned powAlphabetBucketSize = pow(ALPHABETSIZEWITHDOLLAR, bucketSize - 1); //precompute this value, pow is expensive

    //Iterate characterwise and insert them into the buckets
    for (iteration = 1; iteration <= maxLength; ++iteration)
    {

    SEQAN_OMP_PRAGMA(parallel shared(text, bucket, bucketVolume, charCountBuffer, charCountBufferSize, tempBucketInsertCount, bucketInsertCount, pos, qIndex, tempQIndex, buffer))
    {
        String<TValueSize> &threadBktInsertCount = tempBucketInsertCount[omp_get_thread_num()];

        //reset bucketInsertCounter (each thread resets its own counter)
        for (unsigned bucketIndex = 0; bucketIndex < bucketCount; ++bucketIndex)
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

        //for each position in each bucket count the occurences of each character up to this position
        SEQAN_OMP_PRAGMA(for schedule(dynamic))
        for (unsigned bucketIndex = 0; bucketIndex < bucketCount; ++bucketIndex)
        {
            TValueSize const bucketLength = bucketVolume[bucketIndex];
            TValueSize const oldbucketLength = charCountBufferSize[bucketIndex];

            if(oldbucketLength == bucketLength)
            {
                //this bucket did not change, no need to recount
                continue;
            }
            charCountBufferSize[bucketIndex] = bucketLength;

            countCharacters(bucketIndex, charCountBuffer[bucketIndex], bucket[bucketIndex], bucketLength,
                    bucketInsertCount, pos,  ALPHABETSIZE);
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
            TValueSize const startBucketIndex = bucketIndex == 0 ? 0 : bucketInsertCount[bucketIndex - 1];
            TValueSize const endBucketIndex = bucketInsertCount[bucketIndex];
            if (endBucketIndex > startBucketIndex)
                getInsertPositions(startBucketIndex, endBucketIndex, bucketIndex, charCountBuffer, bucket[bucketIndex], charCountBufferSize,
                        qIndex, tempQIndex, pos, threadBktInsertCount, ALPHABETSIZE, powAlphabetBucketSize);
        }

        //accumulate bucket insert count of every thread
        SEQAN_OMP_PRAGMA(for)
        for (unsigned bucketIndex = 0; bucketIndex < bucketCount; ++bucketIndex)
        {
            for (TBktInsertSize t = 1; t < length(tempBucketInsertCount); ++t)
            {
                tempBucketInsertCount[0][bucketIndex] += tempBucketInsertCount[t][bucketIndex];
            }
        }

        //set summarized bucket insert count
        SEQAN_OMP_PRAGMA(single)
        {
            bucketInsertCount[0] = tempBucketInsertCount[0][0];
            for (TValueSize bucketIndex = 1; bucketIndex < bucketCount; ++bucketIndex)
            {
                bucketInsertCount[bucketIndex] = bucketInsertCount[bucketIndex - 1] + tempBucketInsertCount[0][bucketIndex];
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

        //resize Buffer, to avoid resize in parallel section
        SEQAN_OMP_PRAGMA(single nowait)
        {
            for (unsigned bucketIndex = 0; bucketIndex < bucketCount; ++bucketIndex)
            {
                TValueSize const startBucketIndex = bucketIndex == 0 ? 0 : bucketInsertCount[bucketIndex - 1];
                TValueSize const endBucketIndex = bucketInsertCount[bucketIndex];
                resize(buffer[bucketIndex], endBucketIndex - startBucketIndex);
            }
        }

        //Swap current and next qIndex, then sort qIndex, tempQIndex and pos according to the values in qIndex
        SEQAN_OMP_PRAGMA(single nowait)
        {
            swap(qIndex, tempQIndex);
            quickSort(qIndex, tempQIndex, pos, (TValueSize)0, (TValueSize)(length(qIndex) - 1), 0);
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
        for (unsigned bucketIndex = 0; bucketIndex < bucketCount; ++bucketIndex)
        {
            //get new position in ascending order
            TValueSize const startBucketIndex = bucketIndex == 0 ? 0 : bucketInsertCount[bucketIndex - 1];
            TValueSize const endBucketIndex = bucketInsertCount[bucketIndex];
            if (endBucketIndex > startBucketIndex)
                //sorts the buckets section in pos
                doQuickSort(QsortParallel(), pos, sortFunctor, startBucketIndex, endBucketIndex - 1);
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
        for (unsigned bucketIndex = 0; bucketIndex < bucketCount; ++bucketIndex)
        {
            TValueSize const startBucketIndex = bucketIndex == 0 ? 0 : bucketInsertCount[bucketIndex - 1];
            TValueSize const endBucketIndex = bucketInsertCount[bucketIndex];
            if(startBucketIndex < endBucketIndex)
                insertIntoBuckets(startBucketIndex, endBucketIndex, iteration, bucketCount, sentinelSub, text, bucket[bucketIndex], qIndex, tempQIndex, pos, buffer[bucketIndex], bucketVolume[bucketIndex]);
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
    TValueSize index = 0;
    TValueSize sentinelIndex = 0;
    for (unsigned bucketIndex = 0; bucketIndex < bucketCount; ++bucketIndex)
    {
        String<BucketValueType> &currentBucket = bucket[bucketIndex];
        for (TValueSize j = 0; j < length(currentBucket); ++j)
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
inline void preallocateMemory(TText const & text, TBucketCount const & bucketCount, TBucketSize const & bucketSize,
        TAlphabetSize const & alphabetSize, TBucket & bucket, TCharCountBuffer & charCountBuffer)
{
    typedef typename Size<TText>::Type TValueSize;
    TAlphabetSize const alphabetSizeWithDollar = alphabetSize + 1;

#if (SEQAN_ENABLE_DEBUG || SEQAN_ENABLE_TESTING) && SEQAN_ENABLE_PARALLELISM
    double beginAllocTime = omp_get_wtime();
#endif

    String<TValueSize> bucketCapacity;
    resize(bucketCapacity, bucketCount, 0);

    //Compute Bucketsize in parallel
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

//count characters per bucket for each position and store count-values in charCountBuffer
template<typename TBktIndex, typename TCharCountBuffer, typename TBucket, typename TBktLength, typename TBInsertCount,
typename TPos, typename TAlphabetSize>
inline void countCharacters(TBktIndex const bucketIndex, TCharCountBuffer & charCountBuffer, TBucket & currentBucket,
        TBktLength const & bucketLength, TBInsertCount & bucketInsertCount, TPos & pos, TAlphabetSize const & alphabetSize)
{
    typedef typename Value<TBucket>::Type BucketValueType;
    typedef typename Value<TPos>::Type TPosValue;
    typedef typename Value<TPosValue, 2>::Type TValueSize;
    typedef typename Value<TBInsertCount>::Type TBInsertCountValue;

    //Get the first changed index since last iteration. We don't need to recount the previous values
    TBInsertCountValue const startBucketIndex = bucketIndex == 0 ? 0 : bucketInsertCount[bucketIndex - 1];
    TValueSize const & currentPos = pos[startBucketIndex].i2;

    //Initially all counts are 0.
    if (currentPos == 0)
    {
        for (TAlphabetSize j = 0; j < alphabetSize; ++j)
        {
            charCountBuffer[j] = 0;
        }
        BucketValueType const & bktValue = currentBucket[0];
        if (!bktValue.i3)
        {
            ++charCountBuffer[ordValue(bktValue.i2)];
        }
    }


    //iterate over each entry in the bucket, starting at the first changed position
    for (TBktLength index = (currentPos == 0 ? 1 : currentPos); index < bucketLength; ++index)
    {
        TValueSize currentBufferStart = index * alphabetSize;
        TValueSize lastBufferStart = currentBufferStart - alphabetSize;

        typename Iterator<TCharCountBuffer >::Type itCur = iter(charCountBuffer, currentBufferStart);
        typename Iterator<TCharCountBuffer >::Type itLast = iter(charCountBuffer, lastBufferStart);

        //First copy the values from previous entry
        for (TAlphabetSize j = 0; j < alphabetSize; ++j)
        {
            assignValue(itCur, getValue(itLast));
            goNext(itCur);
            goNext(itLast);
            //charCountBuffer[currentBufferStart + j] = charCountBuffer[lastBufferStart + j];
        }

        //Increment count for current character
        BucketValueType const & bktValue = currentBucket[index];
        if (!bktValue.i3)
        {
            ++charCountBuffer[currentBufferStart + ordValue(bktValue.i2)];
        }
    }
}

//For each character inserted into this bucket in the last iteration search the insertposition for the next character of the text
//first count the occurences of the inserted character, which leads to the insertpositions inside the new bucket
//then compute the new bucketindex to insert into
template<typename TBktSize, typename TBktIndex, typename TCharCountBuffer, typename TBucket, typename TCCBSize, typename TQindex, typename TBInsertCount,
typename TPos, typename TAlphabetSize, typename TPowValue>
inline void getInsertPositions(TBktSize const startBucketIndex, TBktSize const endBucketIndex, TBktIndex const bucketIndex, TCharCountBuffer & charCountBuffer, TBucket & currentBucket, TCCBSize & charCountBufferSize,
        TQindex & qIndex, TQindex & tempQIndex, TPos & pos, TBInsertCount & threadBktInsertCount, TAlphabetSize const & alphabetSize, TPowValue const & powAlphabetBucketSize)
{
    typedef typename Value<TBucket>::Type BucketValueType;
    typedef typename Value<TPos>::Type TPosValue;
    typedef typename Value<TPosValue, 2>::Type TValueSize;
    typedef typename Value<TCCBSize>::Type TCCBSizeValue;
    typedef typename Value<TQindex>::Type TQindexValue;

    TAlphabetSize const & alphabetSizeQithDollar = alphabetSize + 1;

    //all positions corresponding to this Bucket
    for (TBktSize i = startBucketIndex; i < endBucketIndex; ++i)
    {
        if (bucketIndex != qIndex[i])
        {
            //this happens only if the text ended in the previous iteration
            //It will be ignored in further iterations
            continue;
        }

        TValueSize const & currentPos = pos[i].i2;
        BucketValueType const & bv = currentBucket[currentPos]; //in previous iteration inserted char
        unsigned const charValue = ordValue(bv.i2);

        //count occurrences in every bucket until we reach current position
        //(For bucketsize>0 it is sufficient to count only the buckets starting with the same character)
        TBktIndex const startBucket = bucketIndex - (bucketIndex % alphabetSizeQithDollar);
        TValueSize charCount = 0;
        for (TBktIndex j = startBucket; j < bucketIndex; ++j)
        {
            TCCBSizeValue const curBufferSize = charCountBufferSize[j];
            if (curBufferSize > 0)
            {
                charCount += charCountBuffer[j][alphabetSize * (curBufferSize - 1) + charValue];
            }
        }
        if (currentPos > 0)
        {
            charCount += charCountBuffer[bucketIndex][alphabetSize * (currentPos - 1) + charValue];
        }

        //remove last char of bucket (shift position (division)) and append new char in front (pow)
        TQindexValue const newBucketIndex = (startBucket / alphabetSizeQithDollar)
                + (powAlphabetBucketSize * (1 + charValue));//1+ to ignore the firstBucket

        pos[i].i2 = charCount;
        tempQIndex[i] = newBucketIndex;
        ++threadBktInsertCount[newBucketIndex];
    }
}

//Insert all previously found characters into this bucket
//first read the character from text, then add the new value to buffer
//after that merge buffer into bucket
template<typename TBktSize, typename TIteration, typename TBktCount, typename TAlphabet, typename TTexts, typename TBucket, typename TQindex,
typename TPos>
inline void insertIntoBuckets(TBktSize const startBucketIndex, TBktSize const endBucketIndex, TIteration const iteration, TBktCount const bucketCount,
        TAlphabet const & sentinelSub, TTexts & text, TBucket & currentBucket,
        TQindex & qIndex, TQindex & tempQIndex, TPos & pos, TBucket & currentBuffer, TBktSize & currentBucketVolume)
{
    typedef typename Value<TBucket>::Type BucketValueType;
    typedef typename Value<TTexts>::Type TText;

    for (TBktSize i = startBucketIndex; i < endBucketIndex; ++i)
    {
       TText const & currentText = text[pos[i].i1];
       bool const isLastCharOfText = iteration == length(currentText);

       //replace with sentinelSub
       TAlphabet const newChar = (isLastCharOfText) ? sentinelSub :
               currentText[length(currentText) - iteration - 1];

       //Insert new values into the buffer. This avoids many move operations within the array
       BucketValueType & cur = currentBuffer[i - startBucketIndex];

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

/*
 * We have 'bufferVolume' new entrys in buffer. Merge them into bucket, so that the sortkeys in value.i1 are ascending.
 */
template<typename BucketValueType, typename TVolume>
void sortBwtBucket(String<BucketValueType> & bucket, TVolume const & bucketVolume,
        String<BucketValueType> const & buffer, TVolume const & bufferVolume)
{
    if (bufferVolume == 0)
        return;

    TVolume insertCount = 0;

    long int curBucketIndex = bucketVolume - 1; //has to be signed value
    TVolume curBufferIndex = bufferVolume - 1;
    TVolume curWriteIndex = bucketVolume + bufferVolume - 1;

    //do a kind of merge sort to insert buffer into bucket, starting with the last values
    while (bufferVolume > insertCount)
    {
        BucketValueType const & curBuf = buffer[curBufferIndex];

        if (curBucketIndex >= 0)
        {
            BucketValueType &curBuk = bucket[curBucketIndex];

            if ((int) (curBuf.i1 - curBufferIndex) <= (int) curBuk.i1) //signed comparison is needed
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
template<typename TAlphabet, typename TFoo1, typename TFoo2, typename TIndex, typename TDepth>
void quickSort(String<TAlphabet> & qIndex, String<TFoo1> & qIndexTemp, String<TFoo2> & pos, TIndex const left,
        TIndex const right, TDepth const depth)
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
