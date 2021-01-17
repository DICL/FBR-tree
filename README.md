# Failure-Atomic Byte-Addressable R-tree for Persistent Memory (FBR-tree)

## Introduction 

In this article, we propose Failure-atomic Byte-addressable R-tree (FBR-tree) that leverages the byte-addressability, persistence, 
and high performance of persistent memory while guaranteeing the crash consistency. 
We carefully control the order of store and cacheline flush instructions and prevent any single store instruction from making 
an FBR-tree inconsistent and unrecoverable. We also develop a non-blocking lock-free range query algorithm for FBR-tree. 
Since FBR-tree allows read transactions to detect and ignore any transient inconsistent states, multiple read transactions can 
concurrently access tree nodes without using shared locks while other write transactions are making changes to them. 
Our performance study shows that FBR-tree successfully reduces the legacy logging overhead and the lock-free range query algorithm 
shows up to 2.6x higher query processing throughput than the shared lock-based crabbing concurrency protocol.
For more details about FBR-tree, please refer to IEEE Transactions on Parallel and Distributed Systems journal
- "[Failure-Atomic Byte-Addressable R-tree for Persistent Memory](https://ieeexplore.ieee.org/abstract/document/9214450)"

## Compilation

```
git clone https://github.com/DICL/FBR-tree
cd FBR-tree
make  
```

## Usage
```
./{exe file} (number_of_INSERT) (number_of_SEARCH) (number_of_Insert_THREADs) (number_of_Search_THREADs) (write_Latency) (k_delta)
```

## Contirubtors
* Soojeong Cho (cristalcho@gmail.com)
* 
