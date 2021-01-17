#include <iostream>
#include <iomanip>
#include <stdint.h>
#define VERSION_MAX UINT8_MAX
#define BIT(N)

template <unsigned int N> // N bytes ( Not Bits )
class MetaData{
  private:
  uint8_t Version_;
  uint8_t Byte[N-1];
  public:
  MetaData(void){
    Version_ = 0;
    for(int i = 0; i < N-1; i++){
      Byte[i] = 0;
    }
  }
//  ~MataData(void){}
  void Init(void){
    Version_ = 0;
    for(int i = 0; i < N-1; i++){
      Byte[i] = 0;
    }
  }
  //Get function
  uint8_t Version(void)const{
    return Version_;
  }
//  uint8_t VersionMax(void)const{
//    return UINT8_MAX;
//  }
  bool IsLeaf(void){
    return Byte[0] & 1;
  }
  bool Bit(int n){
    int Q = n / 8;
    int R = n & 7; // Last 3 bit( 0 ~ 7 )
    //Last bit of Byte[0] is leaf flag.
    if(!(N-2-Q))
      R++;
    return Byte[N-1-1-Q] & (uint8_t)(1 << R);
  }
  
#ifndef INPLACE
  bool SplitCheck(void)const{
    uint8_t cnt = 0;
    for(int i = 0; i < N - 1; i++){
      for(int j = 0; j < 8; j++){
        if(!(i == 0 && j == 0)){
          if(j == 0)
          cnt += ((Byte[i]) & (1));
          else
          cnt += ((Byte[i]>>j) & (1));
        }
      }
    }
    if(cnt >= NODECARD - 1){
      return true;
    }
    return false;
  }
#endif

  bool IsFull(void)const{
    if((Byte[0] & ~(uint8_t)1) != 254){
      return false;
    }
    for (int i = 1; i < N-1; i++){
      if((uint8_t) Byte[i] != 255){
        return false;
      }
    }
    return true;
  }

  /*
  bool IsNearFull(void)const{
    if((Byte[0] & ~(uint8_t)1) != 126){
      return false;
    }
    for (int i = 1; i < N-1; i++){
      if(Byte[i] != 127){
        return false;
      }
    }
    return true;
  }*/

  unsigned int countSetBits(unsigned char n) const{
    unsigned int count = 0;
    while (n) {
        count += n & 1;
        n >>= 1;
    }
    return count;
  }

  bool IsNearFull(void)const{
	int ones = 0;
	char b = Byte[0];

	ones += countSetBits(b>>1); 
    if(ones < 8 - 2)  { // 1-5 
		return false;
    }
    for (int i = 1; i < N-1; i++){
		b = Byte[i];
		ones += countSetBits(b);
		if(ones < ((i+1)*8 - 2) )  {
        	return false;
		}	
	}
    return true;
  }

  //Set function
  void VersionIncr(void){
    Version_++;
  }
  void VersionReset(void){
    Version_ = 0;
  }
  void Leaf(void){
    Byte[0] |= 1;
  }
  void Iter(void){
    Byte[0] &= ~(uint8_t)1;
  }
  void Reset(void){
    Byte[0] &= 1;
    for(int i = 1; i < N-1; i++){
      Byte[i] = 0;
    }
  }
  void Reset(int n){
    int Q = n / 8;
    int R = n & 7; // Last 3 bit( 0 ~ 7 )
//    std::cout << "Q, R: " << N-2-Q << ", " << R << std::endl;
    //Last bit of Byte[0] is leaf flag.
    if(!(N-2-Q))
      R++;
    Byte[N-1-1-Q] &= ~(uint8_t)(1 << R);
  }
  void Set(void){
    Byte[0] &= 1;
    Byte[0] |= 254;
    for(int i = 1; i < N-1; i++){
      Byte[i] = 255;
    }
  }
  void Set(int n){
    int Q = n / 8;
    int R = n & 7; // Last 3 bit( 0 ~ 7 )
    //Last bit of Byte[0] is leaf flag.
    if(!(N-2-Q))
      R++;
    Byte[N-1-1-Q] |= (uint8_t)(1 << R);
  }
  uint8_t* Addr(int nBytes){
    return &Byte[nBytes-1];
  }
  inline int Bit2Byte(int n){
    int Q = n / 8;
    return N - 2 - Q;
  }
  inline int Byte2Atomic(int byte){
    return byte / 8 * 8;
  }
  inline uint8_t* Bit2Addr(int bit){
     return Addr(Byte2Atomic(Bit2Byte(bit)));
  }

  //Utility
  void Print(void)const {
    std::cout << "Version : " << (int)Version_ << std::endl;
    if(Byte[0]&1){
      std::cout << "IsLeaf  : " << "YES" << std::endl;
    } else {
      std::cout << "IsLeaf  : " << "NO " << std::endl;
    }
    for(int i = N-2; i >= 0; i--){
      for(int j = 0; j < 8; j++){
        int P = j + i * 8;
        if(P == 7) break;
        int pos = j + (N-2-i) * 8;
//        std::cout << std::setw(5) << pos << " : ";
        if(P > 7){
          bool var = Byte[i] & (uint8_t)(1 << j);
          std::cout << (int)var << " " ;//std::endl;
        }
        if(P < 7){
          bool var = Byte[i] & (uint8_t)(1 << (j+1));
          std::cout << (int)var << " "; // std::endl;
        }
      }
    }
    std::cout << std::endl; 
  }

};
