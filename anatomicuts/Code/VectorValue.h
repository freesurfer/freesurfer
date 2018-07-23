#ifndef _VectorValue_h_
#define _VectorValue_h_

   template <class T, unsigned int Z=3>
		class VectorValue : public Vector<T,Z>{
			private:
				itk::Vector<T,Z> insideValue;	

			public:
				void SetValue(itk::Vector<T,Z> v )
				{		
					this->insideValue = v;
				}	
				itk::Vector<T,Z> GetValue(){ return this->insideValue; }
		};
   template <class W, class T, unsigned int Z=3>
		class VectorWeighted : public Vector<T,Z>{

      typedef W   ValueType;
      private:
        ValueType insideValue;
			public:
				void SetValue(ValueType v )
				{		
					this->insideValue = v;
				}	
				ValueType GetValue(){ return this->insideValue; }
		};
   template <class T,class V,  unsigned int Z=3>
		class PointValue : public Point<T,Z>{
			private:
				V insideValue;	

			public:
				void SetValue(V v )
				{		
					this->insideValue = v;
				}	
				V GetValue(){ return this->insideValue; }
		};
   template <class W, class T, unsigned int Z=3>
		class PointWeighted : public Point<T,Z>{

      typedef W   ValueType;
      private:
        ValueType insideValue;
			public:
				void SetValue(ValueType v )
				{		
					this->insideValue = v;
				}	
				ValueType GetValue(){ return this->insideValue; }
		};
#endif
