#ifndef Volume3_h
#define Volume3_h


template <typename T>
class Volume3 {
 public:
  Volume3 ( int izX, int izY, int izZ, T iInitialValue );
  ~Volume3 ();

  void SetAll ( T iValue );

  T Get ( int inX, int inY, int inZ );
  void Set ( int inX, int inY, int inZ, T iValue );

  T Get_Unsafe ( int inX, int inY, int inZ );
  void Set_Unsafe ( int inX, int inY, int inZ, T iValue );

 protected:
  T*** mData;
  int mzX;
  int mzY;
  int mzZ;
};


#endif
