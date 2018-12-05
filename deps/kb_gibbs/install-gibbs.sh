#!/bin/bash

curl -O https://raw.githubusercontent.com/arwyer/gibbs_tar/master/gibbs-3.1.tar
tar -xvf gibbs-3.1.tar
cd gibbs-3.1
./configure
sed -i.bak s/-Werror//g Makefile
make
make install
cd ..
