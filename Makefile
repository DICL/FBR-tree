INCLUDES = /path/of/the/hdf5/include/file
LDFLAGS = -Wl,-rpath, /path/of/the/hdf5/lib
LIBS = -lhdf5     
CXX=g++ -std=c++17
CFLAGS=-O3 -w -lrt -lpmemobj -lpthread
all : breakdown breakdown_split_time datasize conc test_hdf5

breakdown:
	$(CXX) -o break-inp-ns1-minlog src/main_for_tests/main.cpp $(CFLAGS) -DINPLACE -DBREAKDOWN 
	$(CXX) -o break-inp-ns2-minlog src/main_for_tests/main.cpp $(CFLAGS) -DINPLACE -DBREAKDOWN -DMULTIMETA -DNS2 
	$(CXX) -o break-inp-ns3-minlog src/main_for_tests/main.cpp $(CFLAGS) -DINPLACE -DBREAKDOWN -DMULTIMETA -DNS3 
	$(CXX) -o break-inp-ns4-minlog src/main_for_tests/main.cpp $(CFLAGS) -DINPLACE -DBREAKDOWN -DMULTIMETA -DNS4 
	$(CXX) -o break-inp-ns1-maxlog src/main_for_tests/main.cpp $(CFLAGS) -DINPLACE -DBREAKDOWN -DFULLLOG 
	$(CXX) -o break-inp-ns2-maxlog src/main_for_tests/main.cpp $(CFLAGS) -DINPLACE -DBREAKDOWN -DMULTIMETA -DNS2 -DFULLLOG 
	$(CXX) -o break-inp-ns3-maxlog src/main_for_tests/main.cpp $(CFLAGS) -DINPLACE -DBREAKDOWN -DMULTIMETA -DNS3 -DFULLLOG 
	$(CXX) -o break-inp-ns4-maxlog src/main_for_tests/main.cpp $(CFLAGS) -DINPLACE -DBREAKDOWN -DMULTIMETA -DNS4 -DFULLLOG 
	$(CXX) -o break-cow-ns1-bitlog src/main_for_tests/main.cpp $(CFLAGS) -DBREAKDOWN -DFULLLOG 
	$(CXX) -o break-cow-ns2-bitlog src/main_for_tests/main.cpp $(CFLAGS) -DBREAKDOWN -DMULTIMETA -DNS2 -DFULLLOG 
	$(CXX) -o break-cow-ns3-bitlog src/main_for_tests/main.cpp $(CFLAGS) -DBREAKDOWN -DMULTIMETA -DNS3 -DFULLLOG 
	$(CXX) -o break-cow-ns4-bitlog src/main_for_tests/main.cpp $(CFLAGS) -DBREAKDOWN -DMULTIMETA -DNS4 -DFULLLOG 

breakdown_split_time:
	$(CXX) -o sbreak-inp-ns1-minlog src/main_for_tests/main.cpp $(CFLAGS) -DINPLACE -DBREAKDOWN -DSPLIT 
	$(CXX) -o sbreak-inp-ns2-minlog src/main_for_tests/main.cpp $(CFLAGS) -DINPLACE -DBREAKDOWN -DMULTIMETA -DNS2 -DSPLIT
	$(CXX) -o sbreak-inp-ns3-minlog src/main_for_tests/main.cpp $(CFLAGS) -DINPLACE -DBREAKDOWN -DMULTIMETA -DNS3 -DSPLIT
	$(CXX) -o sbreak-inp-ns4-minlog src/main_for_tests/main.cpp $(CFLAGS) -DINPLACE -DBREAKDOWN -DMULTIMETA -DNS4 -DSPLIT
	$(CXX) -o sbreak-inp-ns1-maxlog src/main_for_tests/main.cpp $(CFLAGS) -DINPLACE -DBREAKDOWN -DFULLLOG -DSPLIT
	$(CXX) -o sbreak-inp-ns2-maxlog src/main_for_tests/main.cpp $(CFLAGS) -DINPLACE -DBREAKDOWN -DMULTIMETA -DNS2 -DFULLLOG -DSPLIT
	$(CXX) -o sbreak-inp-ns3-maxlog src/main_for_tests/main.cpp $(CFLAGS) -DINPLACE -DBREAKDOWN -DMULTIMETA -DNS3 -DFULLLOG -DSPLIT
	$(CXX) -o sbreak-inp-ns4-maxlog src/main_for_tests/main.cpp $(CFLAGS) -DINPLACE -DBREAKDOWN -DMULTIMETA -DNS4 -DFULLLOG -DSPLIT
	$(CXX) -o sbreak-cow-ns1-bitlog src/main_for_tests/main.cpp $(CFLAGS) -DBREAKDOWN -DFULLLOG -DSPLIT
	$(CXX) -o sbreak-cow-ns2-bitlog src/main_for_tests/main.cpp $(CFLAGS) -DBREAKDOWN -DMULTIMETA -DNS2 -DFULLLOG -DSPLIT
	$(CXX) -o sbreak-cow-ns3-bitlog src/main_for_tests/main.cpp $(CFLAGS) -DBREAKDOWN -DMULTIMETA -DNS3 -DFULLLOG -DSPLIT
	$(CXX) -o sbreak-cow-ns4-bitlog src/main_for_tests/main.cpp $(CFLAGS) -DBREAKDOWN -DMULTIMETA -DNS4 -DFULLLOG -DSPLIT
dclean:
	rm break-*
	rm sbreak-*

datasize:
	$(CXX) -o datasize-cow-ns3-bitlog src/main_for_tests/main.cpp $(CFLAGS) -w -DBREAKDOWN -DMULTIMETA -DNS3 -DFULLLOG
	$(CXX) -o datasize-inp-ns3-minlog src/main_for_tests/main.cpp $(CFLAGS) -w -DINPLACE -DBREAKDOWN -DMULTIMETA -DNS3 

conc:	
	$(CXX) -o conc-inp-ns2-minlog src/main_for_tests/main.cpp $(CFLAGS) -DCONC -INPLACE -DMULTIMETA -DNS2
	$(CXX) -o conc-cow-ns2-bitlog src/main_for_tests/main.cpp $(CFLAGS) -DCONC -DMULTIMETA -DNS2 -DFULLLOG
	$(CXX) -o shared-inp-ns2-minlog src/main_for_tests/main.cpp $(CFLAGS) -DCONC -INPLACE -DSHARED -DMULTIMETA -DNS2
	$(CXX) -o shared-cow-ns2-bitlog src/main_for_tests/main.cpp $(CFLAGS) -DCONC -DSHARED -DMULTIMETA -DNS2 -DFULLLOG
ccclean:
	rm conc?-*
	rm share-*

test_hdf5:  
	$(CXX) -o poissonTest src/main_for_tests/main_poisson.cpp $(INCLUDES) $(LDFLAGS) $(LIBS) $(CFLAGS) -DINPLACE
	$(CXX) -o poissonLockTest src/main_for_tests/main_poisson.cpp $(INCLUDES) $(LDFLAGS) $(LIBS) $(CFLAGS) -DINPLACE -DSHARED  
	$(CXX) -o hdf5test src/main_for_tests/main_hdf5.cpp $(INCLUDES) $(LDFLAGS) $(LIBS) $(CFLAGS) -DINPLACE -DCONC 
	$(CXX) -o Lockhdf5test src/main_for_tests/main_hdf5.cpp $(INCLUDES) $(LDFLAGS) $(LIBS) $(CFLAGS) -DINPLACE -DCONC -DSHARED 
	$(CXX) -o onlySearch src/main_for_tests/main_hdf5.cpp $(INCLUDES) $(LDFLAGS) $(LIBS) $(CFLAGS) -DINPLACE -DCONC -DOS 

hdclean:
	rm poisson*
	rm hdf5test
	rm onlySearch
	rm Lockhdf5test

#taxi : 83409385, HDF5 : 4486 * 100
