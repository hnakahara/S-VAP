#!/bin/bash

# error of libcrypto.so.1.1
for i in `find . -name libcrypto.so.1.1`
do
  echo $i
  replicant=${i%.1}.0.0
  echo $replicant
  if [ -e $i ] && [ ! -e $replicant ]; then
    echo "copy execution"
    cp -n $i $replicant
  else
    echo "file already exist"
  fi
done
