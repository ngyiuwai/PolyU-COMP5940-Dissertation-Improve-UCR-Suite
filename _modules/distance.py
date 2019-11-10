"""
This module is for calculating distance/ lower bound distance between two sequences with equal length.
Functions:  squareEuclidean (normQuery, seq, seqMean, seqStd, order, bestSoFar) --> float
            squareDTW (query, rawSeq, rawSeqMean, rawSeqStd, order, bestSoFar, constrain) --> float
            LB_Kim() -->
            LB_Keogh() -->
            LB_Keogh2() -->
"""


def squareEuclidean(query: list, seq: list, seqMean: float, seqStd: float, order: list, bestSoFar: float):
    """
    This function is for calculating the square of Euclidean Distance between two sequences.
    If the square distance is greater than "best-so-far" square distance, we should assume that distance is infinity.
    IMPORTANT:  "Early Abandoning Z-Normalization" technique of UCR Suite is applied here.
                "Reordering Early Abandoning" technique of UCR Suite is applied here. We assume input sequence and query are reordered.
    ================================================================
    Input:
      query:    z-normalized query.
      seq:      Non-normalized subsequence.
      seqMean:  Mean of sequence, which is for z-normalization.
      seqStd:   Standard derivation of sequence, which is for z-normalization.
      order:    The descending order of query's absolute z-normalized value.
      bestSoFar:  stop calculating square ED if the output is greater than best-so-far. Best-so-far is a squared ED too.
    Output:
      squareDistance: float
    """
    distance = 0
    for i in order:
        d = (seq[i] - seqMean) / seqStd - query[i]
        distance = distance + d*d
        if (distance > bestSoFar):
            distance = float('inf')
            break
    return distance


def dynamicTimeWraping(query: list, seq: list, seqMean: float, seqStd: float, cumLB: list, scBand: int, bestSoFar: float):
    """
    This function is for calculating the Dynamic Time Wrapping square distance between two sequences.
    If the square distance is greater than "best-so-far" square distance, we should assume that distance is infinity.
    IMPORTANT:  "Early Abandoning Z-Normalization" technique of UCR Suite is applied here.
                We cannot apply "Reordering Early Abandoning" technique because of the nature of DTW. 
                Such technique was applied in Lower Bound. But still we can use "Sequential early abandoning" using cumLB[].
    ================================================================
    Input:
      query:    z-normalized query.
      seq:      Non-normalized subsequence.
      seqMean:  Mean of sequence, which is for z-normalization.
      seqStd:   Standard derivation of sequence, which is for z-normalization.
      cumLB:    Lower bound distance of each points calculated in LB_Keogh
      scBand:   Sakoe-Chiba warpping band
      bestSoFar:  stop calculating square ED if the output is greater than best-so-far. Best-so-far is a squared ED too.
    Output:
      squareDistance: float
    ================================================================
    ++++++ It is an approximation of DTW. ++++++
    ++++++ The aim of this project is to find a lower bound of LB_Keogh. Speed up the actual DTW process is not our concern. ++++++
    ++++++ Please refer to function "dynamicTimeWraping_true()" for real DTW, which is very very slow. ++++++
    """
    route = [len(query) - 1, len(query) - 1]
    # DWT is a path searching algorthm. Our goal is to move the point from end to begin
    # x-coordinate = position in query
    # y-coordinate = position in sequence
    # The distance at ending point is always included.
    d = (seq[route[1]] - seqMean) / seqStd - query[route[0]]
    distance = d*d
    lb = cumLB[route[0]]
    while (route[0] > 0) or (route[1] > 0):
        lb = lb + cumLB[route[0]]
        if (lb > bestSoFar) or (distance > bestSoFar):
            distance = float('inf')
            break
        # Search distance from the 3 neighbor of route = (i, j)
        # i.e.  distance(i-1, j)
        #       distance(i, j-1)
        #       distance(i-1, j-1)
        # Reminder: i = position in query, j = position in sequence
        d1 = abs((seq[route[1]] - seqMean) / seqStd - query[route[0] - 1])
        d2 = abs((seq[route[1] - 1] - seqMean) / seqStd - query[route[0]])
        d3 = abs((seq[route[1] - 1] - seqMean) / seqStd - query[route[0] - 1])

        # Start path searching
        if (route[0] == 0):
            distance = distance + d2*d2
            route[1] = route[1] - 1
        elif (route[1] == 0):
            distance = distance + d1*d1
            route[0] = route[0] - 1

        elif (route[0] - route[1]) >= scBand:
            if (d3 < d1):
                distance = distance + d3*d3
                route[0] = route[0] - 1
                route[1] = route[1] - 1
            else:
                distance = distance + d1*d1
                route[0] = route[0] - 1
        elif (route[1] - route[0]) >= scBand:
            if (d3 < d2):
                distance = distance + d3*d3
                route[0] = route[0] - 1
                route[1] = route[1] - 1
            else:
                distance = distance + d2*d2
                route[1] = route[1] - 1

        elif (d3 < d1) and (d3 < d2):
            distance = distance + d3*d3
            route[0] = route[0] - 1
            route[1] = route[1] - 1
            pass
        elif (d2 < d1):
            distance = distance + d2*d2
            route[1] = route[1] - 1
            pass
        else:
            distance = distance + d1*d1
            route[0] = route[0] - 1
            pass

    return distance


def dynamicTimeWraping_true(query: list, seq: list, seqMean: float, seqStd: float, cumLB: list, scBand: int, bestSoFar: float):
    """
    Real DTW which use dynamic programming. Very very slow for long query.
    Python is not designed of dynamic programming. Use C/C++ instead...
    ================================================================
    """

    # Distance at current position, i.e. ending point, is always included.
    posQ = len(query) - 1
    posS = len(seq) - 1
    d = (seq[posS] - seqMean) / seqStd - query[posQ]
    distCurrent = d*d

    lb = sum(cumLB[0:posQ]) + distCurrent

    if lb > bestSoFar:
        return float('inf')

    # Return to uppper level if current postion is (0, 0)
    if (posQ == 0) and (posS == 0):
        return distCurrent
    else:
        # If current position hits the left bound, next position should be bottom.
        if (posS == 0):
            query.pop()
            distBottom = dynamicTimeWraping_true(
                query, seq, seqMean, seqStd, cumLB, scBand, bestSoFar)
            return distBottom + distCurrent
        # If current position hits the bottom bound, next position should be left.
        elif (posQ == 0):
            seq.pop()
            distLeft = dynamicTimeWraping_true(
                query, seq, seqMean, seqStd, cumLB, scBand, bestSoFar)
            return distLeft + distCurrent

        # If current position hits the left SC band, next postion should be diagonal or bottom
        if (posQ - (posS - 1)) > scBand:
            qryDiag = list(query)
            qryDiag.pop()
            seqDiag = list(seq)
            seqDiag.pop()
            qryBottom = list(query)
            qryBottom.pop()
            seqBottom = list(seq)
            distDiag = dynamicTimeWraping_true(
                qryDiag, seqDiag, seqMean, seqStd, cumLB, scBand, bestSoFar)
            distBottom = dynamicTimeWraping_true(
                qryBottom, seqBottom, seqMean, seqStd, cumLB, scBand, bestSoFar)
            if (distBottom < distDiag):
                query.pop()
                return distBottom + distCurrent
            else:
                query.pop()
                seq.pop()
                return distDiag + distCurrent

        # If current position hits the bottom SC band, next postion should be diagonal or left
        elif (posS - (posQ - 1)) > scBand:
            qryLeft = list(query)
            seqLeft = list(seq)
            seqLeft.pop()
            qryDiag = list(query)
            qryDiag.pop()
            seqDiag = list(seq)
            seqDiag.pop()
            distLeft = dynamicTimeWraping_true(
                qryLeft, seqLeft, seqMean, seqStd, cumLB, scBand, bestSoFar)
            distDiag = dynamicTimeWraping_true(
                qryDiag, seqDiag, seqMean, seqStd, cumLB, scBand, bestSoFar)
            if (distLeft < distDiag):
                seq.pop()
                return distLeft + distCurrent
            else:
                query.pop()
                seq.pop()
                return distDiag + distCurrent

        # If current position dose not hit any bound, next position could be diagonal, left or bottom
        else:
            qryLeft = list(query)
            seqLeft = list(seq)
            seqLeft.pop()
            qryDiag = list(query)
            qryDiag.pop()
            seqDiag = list(seq)
            seqDiag.pop()
            qryBottom = list(query)
            qryBottom.pop()
            seqBottom = list(seq)
            distLeft = dynamicTimeWraping_true(
                qryLeft, seqLeft, seqMean, seqStd, cumLB, scBand, bestSoFar)
            distDiag = dynamicTimeWraping_true(
                qryDiag, seqDiag, seqMean, seqStd, cumLB, scBand, bestSoFar)
            distBottom = dynamicTimeWraping_true(
                qryBottom, seqBottom, seqMean, seqStd, cumLB, scBand, bestSoFar)
            if (distLeft < distDiag) and (distLeft < distBottom):
                seq.pop()
                return distLeft + distCurrent
            elif (distDiag < distBottom):
                query.pop()
                seq.pop()
                return distDiag + distCurrent
            else:
                query.pop()
                return distBottom + distCurrent


def LB_Kim(query: list, seq: list, seqMean: float, seqStd: float, cumLB: list, bestSoFar: float):
    """
    This function is for calculating the LB_Kim between two sequences.
    If the square distance is greater than "best-so-far" square distance, we should assume that distance is infinity.
    LB_Kim says the front, back, top, bottem are lower bounds.
    But the original UCR Suite (written in cpp) suggest than z-normalization the top and bottom cannot give siginifant benefits.
    Hence, we will compute distance at 3 points at front and 3 points at back only.
    IMPORTANT:  "Early Abandoning Z-Normalization" technique of UCR Suite is applied here.
                "Reordering Early Abandoning" technique of UCR Suite is applied here. We assume input sequence and query are reordered.
    ================================================================
    Input:
      query:    z-normalized query.
      seq:      Non-normalized subsequence.
      seqMean:  Mean of sequence, which is for z-normalization.
      seqStd:   Standard derivation of sequence, which is for z-normalization.
      cumLB:    Lower bound distance of each points calculated. Initially empty
      bestSoFar:  stop calculating square ED if the output is greater than best-so-far. Best-so-far is a squared ED too.
    Output:
      squareDistance: float
      cumLB:    Lower bound distance of each points calculated in LB_Kim.
    """

    # 1st point at front and back
    sfront0 = (seq[0] - seqMean) / seqStd
    dfront0 = sfront0 - query[0]
    dfront0 = dfront0 * dfront0

    sback0 = (seq[len(seq)-1] - seqMean) / seqStd
    dback0 = sback0 - query[len(query)-1]
    dback0 = dback0 * dback0

    cumLB[0] = dfront0
    cumLB[len(query)-1] = dback0
    distance = dfront0 + dback0

    if (distance > bestSoFar):
        distance = float('inf')
        return distance, cumLB

    # 2nd points at front and back
    sfront1 = (seq[1] - seqMean) / seqStd
    dfront1 = min(abs(sfront1 - query[1]),
                  abs(sfront1 - query[0]))
    dfront1 = dfront1 * dfront1
    dfront1 = min(dfront1, dfront0)

    sback1 = (seq[len(query)-2] - seqMean) / seqStd
    dback1 = min(abs(sback1 - query[len(query)-1]),
                 abs(sback1 - query[len(query)-2]))
    dback1 = dback1 * dback1
    dback1 = min(dback1, dback0)

    cumLB[1] = dfront1
    cumLB[len(query)-2] = dback1
    distance = distance + dfront1 + dback1
    if (distance > bestSoFar):
        distance = float('inf')
        return distance, cumLB

    # 3rd points at front and back
    sfront2 = (seq[2] - seqMean) / seqStd
    dfront2 = min(abs(sfront2 - query[2]),
                  abs(sfront2 - query[1]))
    dfront2 = dfront2 * dfront2
    dfront2 = min(dfront2, dfront1)

    sback2 = (seq[len(query)-3] - seqMean) / seqStd
    dback2 = min(abs(sback2 - query[len(query)-1]),
                 abs(sback2 - query[len(query)-2]))
    dback2 = dback2 * dback2
    dback2 = min(dback2, dback1)

    cumLB[2] = dfront2
    cumLB[len(query)-3] = dback2
    distance = distance + dfront2 + dback2
    if (distance > bestSoFar):
        distance = float('inf')
        return distance, cumLB

    return distance, cumLB


def LB_Keogh(queryUp: list, queryLow: list, seq: list, seqMean: float, seqStd: float, order: list, cumLB: list, bestSoFar: float):

    distance = 0
    for i in order:
        d = 0
        normSeq = (seq[i] - seqMean) / seqStd
        if (queryUp[i] < normSeq):
            d = normSeq - queryUp[i]
            distance = distance + d*d
        if (queryLow[i] > normSeq):
            d = normSeq - queryLow[i]
            distance = distance + d*d
        cumLB[i] = d*d

        if (distance > bestSoFar):
            distance = float('inf')
            break

    return distance, cumLB


def LB_Keogh2(query: list, seqUp: list, seqLow: list, seqMean: float, seqStd: float, order: list, cumLB: list, bestSoFar: float):

    distance = 0
    for i in order:
        d = 0
        normSeqUp = (seqUp[i] - seqMean) / seqStd
        normSeqLow = (seqLow[i] - seqMean) / seqStd
        if (query[i] > normSeqUp):
            d = normSeqUp - query[i]
            distance = distance + d*d
        if (query[i] < normSeqLow):
            d = normSeqLow - query[i]
            distance = distance + d*d
        cumLB[i] = max(cumLB[i], d*d)

        if (distance > bestSoFar):
            distance = float('inf')
            break

    return distance, cumLB
