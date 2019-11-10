from _modules import distance
from _modules import sort
import os
import time
import json
import matplotlib.pyplot as plt

# Read Query

script_dir = os.path.dirname(__file__)
read_path = os.path.join(script_dir, 'rawData/query.json')
textfile = open(read_path, "r")
query = json.loads(textfile.read())
textfile.close()

# Normalize query
queryMean = sum(query) / len(query)
queryStd = 0
for i in range(0, len(query)):
    queryStd = queryStd + query[i]*query[i]
queryStd = (queryStd/len(query) - queryMean*queryMean) ** 0.5

normQuery = []
cumLB = []  # cumulative lower bound, which is not used in raw DTW
for i in range(0, len(query)):
    normQuery.append((query[i] - queryMean) / queryStd)
    cumLB.append(0)

# Define best-so-far, as we are returning query distance < bsf
# It is k-nearest neighor seacrch, so bsf is initially set to INF.
bsf = float("inf")
k = 5
scBand = 5
print('query mean   =', queryMean)
print('query std    =', queryStd)
print('query len    =', len(query))
print('best-so-far  =', bsf)
print('kNN neighbor =', k)
_, sortingOrder = sort.bubbleSort(normQuery)
print('Sakoe-Chiba  =', scBand)
print('ordering     =', sortingOrder)

plt.plot(normQuery)
plt.title('Normalized Query')
plt.show()

### Step 1: Read all raw data ###
# For fair comparision, ALL data are loaded in memory.

print('......', time.ctime(), ' start loading data ......')

script_dir = os.path.dirname(__file__)
read_path = os.path.join(script_dir, 'rawData/data.json')
textfile = open(read_path, "r")
data = json.loads(textfile.read())
textfile.close()
print('......', time.ctime(), ' finsih loading data ......')

# Step 2: Run ED

input('Please Enter to run similarity search by UCR_DTW')
print('......', time.ctime(), ' start running UCR_DTW ......')

match = []

# First k subsequence are assume to be nearest neigbour.

# Prepare sum of x and sum of x^2 for fast mean and std calculation.
# 1st subsequence : Full computation
subseq = data[0:len(query)]
subSeqSum = 0
subSeqSum2 = 0
for i in range(0, len(query)):
    subSeqSum = subSeqSum + data[i]
    subSeqSum2 = subSeqSum2 + data[i]*data[i]
subseqMean = subSeqSum / len(query)
subseqStd = (subSeqSum2 / len(query) - subseqMean * subseqMean) ** 0.5
q = list(normQuery)
s = list(subseq)
dist = distance.dynamicTimeWraping(
    q, s, subseqMean, subseqStd, cumLB, scBand, bsf)
match.append([dist, i])

# 2nd to len(query) subsequence : Online update
for i in range(1, k):
    old = subseq.pop(0)
    new = data[i+len(query)]
    subseq.append(new)
    subSeqSum = subSeqSum - old + new
    subSeqSum2 = subSeqSum2 - old*old + new*new
    subseqMean = subSeqSum / len(query)
    subseqStd = (subSeqSum2 / len(query) - subseqMean * subseqMean) ** 0.5
    q = list(normQuery)
    s = list(subseq)
    dist = distance.dynamicTimeWraping(
        q, s, subseqMean, subseqStd, cumLB, scBand, bsf)
    match.append([dist, i])

match.sort()

# Replace k-nearest neigbour if dist is shorter.
bsf = match[k-1][0]

for i in range(k, len(data)-len(query)):
    old = subseq.pop(0)
    new = data[i+len(query)]
    subseq.append(new)
    subSeqSum = subSeqSum - old + new
    subSeqSum2 = subSeqSum2 - old*old + new*new
    subseqMean = subSeqSum / len(query)
    subseqStd = (subSeqSum2 / len(query) - subseqMean * subseqMean) ** 0.5
    q = list(normQuery)
    s = list(subseq)
    dist = distance.dynamicTimeWraping(
        q, s, subseqMean, subseqStd, cumLB, scBand, bsf)
    if (dist < float('inf')):
        replaced = match.pop()
        match.append([dist, i])
        match.sort()
        bsf = match[k-1][0]
        # print('Replacing #', replaced[1], ' by #', i)

print('......', time.ctime(), ' finish running UCR_DTW ......')
print('kNN: ')
print(match)
# Plot results


for i in range(0, k):
    plt.plot(data[match[i][1]:match[i][1]+len(query)])
    plt.title(match[i][1])
    plt.show()
