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

#ifndef CORE_INCLUDE_SEQAN_INDEX_QUICKSORT_H_
#define CORE_INCLUDE_SEQAN_INDEX_QUICKSORT_H_

namespace seqan
{

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================
struct qsortParallel_;
typedef Tag<qsortParallel_> QsortParallel;

struct qsortSequential_;
typedef Tag<qsortSequential_> QsortSequential;
// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// sort SA according to sortKeys
// choose between quickSort and insertionSort
//

/**
.Function.doQuickSort:
..summary:sorts the input
..description:
The Input gets sortet via quicksort, according to the sortFunctor, if the range between start and end gets too small, insertionSort is used
..doQuickSort(tag, toSort, sortFunctor, start, end)
..param.tag: Defines whether quicksort uses openMp Tasks or sequential computation
..param.toSort: the input which gets sorted
..param.sortFunctor: Functor returns the sortValue for a given index
..param.start: startIndex to sort
..param.end: endIndex to sort
*/
template<typename TToSort, typename TSortFunctor, typename TIndex>
void doQuickSort(QsortSequential const tag, TToSort & toSort, TSortFunctor & sortFunctor, TIndex const start, TIndex const end)
{
    if (end - start > 15)
        quickSort(tag, toSort, sortFunctor, start, end, 0);
    else
        insertionSort(toSort, sortFunctor, start, end);
}

/**
.Function.doQuickSort:
..summary:sorts the input
..description:
The Input gets sortet via quicksort, according to the sortFunctor, if the range between start and end gets too small, insertionSort is used
..doQuickSort(tag, toSort, sortFunctor, start, end)
..param.tag: Defines whether quicksort uses openMp Tasks or sequential computation
..param.toSort: the input which gets sorted
..param.sortFunctor: Functor returns the sortValue for a given index
..param.start: startIndex to sort
..param.end: endIndex to sort
*/
template<typename TToSort, typename TSortFunctor, typename TIndex>
void doQuickSort(QsortParallel const tag, TToSort & toSort, TSortFunctor sortFunctor, TIndex const start, TIndex const end)
{
    if (end - start > 15)
        quickSort(tag, toSort, sortFunctor, start, end, 0);
    else
        insertionSort(toSort, sortFunctor, start, end);
}

/**
.Function.doQuickSort:
..summary:sorts the input
..description:
The Input gets sortet via quicksort, according to the sortFunctor, if the range between start and end gets too small, insertionSort is used
..doQuickSort(tag, toSort, sortFunctor, start, end, depth)
..param.tag: Defines whether quicksort uses openMp Tasks or sequential computation
..param.toSort: the input which gets sorted
..param.sortFunctor: Functor returns the sortValue for a given index
..param.start: startIndex to sort
..param.end: endIndex to sort
..param.depth: depth of recursion, should be 0 initially
*/
template<typename TToSort, typename TSortFunctor, typename TIndex>
void doQuickSort(QsortParallel const tag, TToSort & toSort, TSortFunctor sortFunctor, TIndex const start, TIndex const end,
       long const depth)
{
    if (end - start > 15)
        quickSort(tag, toSort, sortFunctor, start, end, depth);
    else
        insertionSort(toSort, sortFunctor, start, end);
}

template<typename TToSort, typename TSortFunctor, typename TIndex>
void doQuickSort(QsortSequential const tag, TToSort & toSort, TSortFunctor sortFunctor, TIndex const start, TIndex const end,
       long const depth)
{
    if (end - start > 15)
        quickSort(tag, toSort, sortFunctor, start, end, depth);
    else
        insertionSort(toSort, sortFunctor, start, end);
}

//gets the median of the three values
//returns 0, 1 or 2 for a, b or c
template<typename Val>
int medianOfThree(Val const a, Val const b, Val const c)
{
    if (a == b)
        return 0;
    if (a == c || b == c)
        return 2;
    if (a < b)
    {
        if (b < c)
            return 1;
        else if (a < c)
            return 2;
        else
            return 0;
    }
    else if (b > c)
        return 1;
    else if (a < c)
        return 0;
    else
        return 2;
}

//finds and returns the pivot element within the given range using medianOfThree
template<typename TResult, typename TToSort, typename TSortFunctor, typename TIndex>
TResult _findPivot(TToSort & toSort, TSortFunctor sortFunctor, TIndex const left, TIndex const right)
{
    typedef typename TSortFunctor::result_type TSortKey;

    TIndex const halfBucketSize = (right - left) / 2;
    TSortKey pivot = sortFunctor(right);
    TSortKey const pivotB = sortFunctor(left);
    TSortKey const pivotC = sortFunctor(left + halfBucketSize);
    int const medianNumber = medianOfThree(pivot, pivotB, pivotC);

    if (medianNumber != 0)
    {
        pivot = medianNumber == 1 ? pivotB : pivotC;
        TIndex const swapIndex = medianNumber == 1 ? left : left + halfBucketSize;
        std::swap(toSort[swapIndex], toSort[right]);
    }
    return pivot;
}

//partitions the elements according the pivot element
template<typename TResult, typename TPivot, typename TToSort, typename TSortFunctor, typename TIndex>
TResult _partition(TToSort & toSort, TPivot const pivot, TSortFunctor sortFunctor, TIndex const left, TIndex const right)
{
    //lomutos partitioning scheme
    TIndex leftTmp = left;
    TIndex rightTmp = left;
    while (rightTmp < right)
    {
        if (sortFunctor(rightTmp) <= pivot)
        {
            std::swap(toSort[leftTmp], toSort[rightTmp]);

            leftTmp++;
        }
        rightTmp++;
    }

    std::swap(toSort[leftTmp], toSort[rightTmp]);

    return leftTmp;
}

//Sort the input with Quicksort. Partition with the pivot element, then do a recursion
template<typename TToSort, typename TSortFunctor, typename TIndex>
void quickSort(QsortSequential const tag, TToSort & toSort, TSortFunctor sortFunctor, TIndex const left, TIndex const right, long depth)
{
    if(depth > 1000)
    {
        insertionSort(toSort, sortFunctor, left, right);
        return;
    }

    typedef typename TSortFunctor::result_type TSortKey;

    TSortKey const pivot = _findPivot<TSortKey>(toSort, sortFunctor, left, right);

    TIndex leftTmp = _partition<TIndex>(toSort, pivot, sortFunctor, left, right);

    if (right > leftTmp + 1)
    {
        doQuickSort(tag, toSort, sortFunctor, leftTmp + 1, right, depth + 1);
    }
    do
    {
        --leftTmp;
    }
    while (leftTmp > left && sortFunctor(leftTmp) == pivot);
    if (leftTmp > left)
    {
        doQuickSort(tag, toSort, sortFunctor, left, leftTmp, depth + 1);
    }
}

//Sort the input with Quicksort. Partition with the pivot element, then do a recursion
template<typename TToSort, typename TSortFunctor, typename TIndex>
void quickSort(QsortParallel const tag, TToSort & toSort, TSortFunctor sortFunctor, TIndex const left, TIndex const right,
        long depth)
{
    if (depth > 10)
    {
        //Avoid the creation of too many tasks. Their overhead would make it very slow
        quickSort(QsortSequential(), toSort, sortFunctor, left, right, depth);
        return;
    }

    typedef typename TSortFunctor::result_type TSortKey;

    TSortKey const pivot = _findPivot<TSortKey>(toSort, sortFunctor, left, right);

    TIndex leftTmp = _partition<TIndex>(toSort, pivot, sortFunctor, left, right);

    SEQAN_OMP_PRAGMA(task firstprivate(leftTmp, right, depth, sortFunctor) shared(toSort))
    {
        if (right > leftTmp + 1)
        {
            doQuickSort(tag, toSort, sortFunctor, leftTmp + 1, right, depth + 1);
        }
    }
    do
    {
        --leftTmp;
    }
    while (leftTmp > left && sortFunctor(leftTmp) == pivot);
    if (leftTmp > left)
    {
        doQuickSort(tag, toSort, sortFunctor, left, leftTmp, depth + 1);
    }
}

//sort the input with insertionSort
template<typename TToSort, typename TSortFunctor, typename TIndex>
void insertionSort(TToSort & toSort, TSortFunctor sortFunctor, TIndex const left, TIndex const right)
{
    typedef typename Value<TToSort>::Type TToSortVal;
    typedef typename TSortFunctor::result_type TSortKey;

    TIndex j = 0;
    TSortKey sortkey;
    TToSortVal tmp;

    for (TIndex i = left + 1; i <= right; ++i)
    {
        sortkey = sortFunctor(i);
        tmp = toSort[i];
        j = i - 1;

        int c = 1;
        while ((j >= left) && (sortFunctor(j) > sortkey))
        {
            toSort[j + 1] = toSort[j];
            if (j == 0)
            {
                c = 0;
                break; //underflow with unsigned
            }
            --j;
        }
        toSort[j + c] = tmp;
    }
}

}  // namespace seqan

#endif  // #ifndef CORE_INCLUDE_SEQAN_INDEX_QUICKSORT_H_
