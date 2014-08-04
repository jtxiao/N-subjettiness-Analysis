N-subjettiness-Analysis
=======================

Instruction:

* Need to install fastjet (for instance 3.0.6)

curl -O http://fastjet.fr/repo/fastjet-3.0.6.tar.gz 

tar zxvf fastjet-3.0.6.tar.gz

cd fastjet-3.0.6/

./configure --prefix=$PWD/../fastjet-install

make -j 20

make check  

make install

cd ..


* Fastjet contrib
 
wget http://fastjet.hepforge.org/contrib/downloads/fjcontrib-1.013.tar.gz

tar xvfz fjcontrib-1.013.tar.gz

cd fjcontrib-1.013/

./configure --fastjet-config=$PWD/../fastjet-3.0.6/fastjet-config                  

make -j 20

make check  
  
make install

* HiForestAnalysis helper class

git clone git@github.com:CmsHI/HiForestAnalysis.git

