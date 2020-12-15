#!/bin/bash

#cat foldchnge_fragment_norm* > foldchnge_6_norm.txt

cat foldchnge_fragment_noimp.txt foldchnge_fragment_0.txt foldchnge_fragment_minDet.txt foldchnge_fragment_minProb.txt foldchnge_fragment_svd.txt foldchnge_fragment_bpca.txt \
foldchnge_fragment_knn.txt foldchnge_fragment_lls.txt foldchnge_fragment_rf.txt > foldchnge_fragment_9_norm.txt
