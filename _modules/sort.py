"""
This module is for sorting and reordering sequence.
Functions:  sort (query) --> (sortedQuery, sortingOrder)
            reorder (sequence, sortingOrder) --> sortedSequence
"""

def bubbleSort (seq: list):
    """
    This function is for sorting a z-normalized sequence descending by absolute z-normalization value.
    We need to record how we reorder query, so sequences could be reordered in the same way for comparision.
    For convenient, bubble sort is used.
    ================================================================
    Input:
      seq:  z-normalized subsequence
    Output:
      (sortedSeq, sortingOrder): (list, list)
    """
    sortedSeq = list(seq)
    sortingOrder = list(range(len(seq)))

    for i in range(0, len(seq)):
      for j in range(0, len(seq) - i - 1):
        if abs(sortedSeq[j]) < abs(sortedSeq[j+1]):
          sortedSeq[j], sortedSeq[j+1] = sortedSeq[j+1], sortedSeq[j]
          sortingOrder[j], sortingOrder[j+1] = sortingOrder[j+1], sortingOrder[j]

    return sortedSeq, sortingOrder

      