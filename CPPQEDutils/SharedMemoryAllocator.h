// Copyright András Vukics 2006–2022. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#ifndef CPPQEDCORE_UTILS_SHAREDMEMORYALLOCATOR_H_INCLUDED
#define CPPQEDCORE_UTILS_SHAREDMEMORYALLOCATOR_H_INCLUDED


namespace cppqedutils {

  
template<typename> class SharedMemoryAllocator;

  
template <>
class SharedMemoryAllocator<void>
{
public:
  typedef void* pointer;
  typedef const void* const_pointer;
  typedef void value_type;
  template<class U> 
  struct rebind {typedef SharedMemoryAllocator<U> other;};
};


template<typename T>
class SharedMemoryAllocator
{
public:
  // type definitions
  typedef       T   value_type;
  typedef       T*  pointer;
  typedef const T*  const_pointer;
  typedef       T&  reference;
  typedef const T&  const_reference;
  typedef size_t    size_type;
  typedef ptrdiff_t difference_type;
    
  // rebind SharedMemoryAllocator to type U
  template<class U>
  struct rebind {typedef SharedMemoryAllocator<U> other;};

  // return address of values
  pointer address(reference value) const {return &value;}
  const_pointer address(const_reference value) const {return &value;}
    
  // constructors and destructor
  SharedMemoryAllocator() throw() : ptrToData_(0) {std::cerr<<"Constructor: ptrToData_="<<ptrToData_<<std::endl;}
  SharedMemoryAllocator(const pointer x, size_t size) : ptrToData_(x), size_(size) {std::cerr<<"Constructor: ptrToData_="<<ptrToData_<<std::endl;}
  SharedMemoryAllocator(const SharedMemoryAllocator& a) throw() : ptrToData_(a.ptrToData_) {};
  template <class U> SharedMemoryAllocator(const SharedMemoryAllocator<U>&) throw() {}
  ~SharedMemoryAllocator() throw() {}
    
  // return maximum number of elements that can be allocated
  size_type max_size () const throw() {return std::numeric_limits<size_t>::max() / sizeof(T);}
    
  // allocate but do not initialize num elements of type T
  pointer allocate (size_type num, const void* =0)
  {
    std::cerr<<"Allocation: ptrToData_="<<ptrToData_<<std::endl;
    if (ptrToData_==0)
      {
        // print message and allocate memory with global new
        pointer ret=(pointer)(::operator new(num*sizeof(T)));
        std::cerr<<"allocate "<<num<<" element(s) of size "<<sizeof(T)<<" allocated at: "<<(void*)ret<<std::endl;
        return ret;
      }
    else
      {
        std::cerr << "Allocation: data passed"<<std::endl;
        if (num!=size_) ERROR!!!
        pointer ret=ptrToData_;
      }
  }
    
  // initialize elements of allocated storage p with value value
  void construct (pointer p, const T& value)
  { 
    std::cerr<<"Initialization: ptrToData_="<<ptrToData_<<std::endl;
    if (ptrToData_==0)
      {
        // initialize memory with placement new
        std::cerr << "Initialization"<<std::endl;
        new((void*)p)T(value);
      }
    else
      {
        std::cerr << "Initialization: data passed"<<std::endl;
        p=ptrToData_; //in fact not necessary
      }
  }
    
  // destroy elements of initialized storage p
  void destroy (pointer p)
  {
    std::cerr<<"Destroy: ptrToData_="<<ptrToData_<<std::endl;
    if (ptrToData_==0)
      {
        // destroy objects by calling their destructor
        p->~T();
      }
    else
      {
        std::cerr << "Destroy: data passed"<<std::endl;
      }
  }
    
  // deallocate storage p of deleted elements
  void deallocate (pointer p, size_type num)
  {
    std::cerr<<"Deallocate: ptrToData_="<<ptrToData_<<std::endl;
    if (ptrToData_==0)
      {
        // print message and deallocate memory with global delete
        std::cerr<<"deallocate "<<num<<" element(s) of size "<<sizeof(T)<<" at: "<<(void*)p<<std::endl;
        ::operator delete((void*)p);
      }
    else
      {
        std::cerr << "Deallocate: data passed"<<std::endl;
        p=0;
      }
  }

private:
  pointer const ptrToData_;
  const size_t  size_;
     
}; //SharedMemoryAllocator
  
// return that all specializations of this SharedMemoryAllocator are interchangeable
template <class T1, class T2>
bool operator== (const SharedMemoryAllocator<T1>&, const SharedMemoryAllocator<T2>&) throw() {return true;}
template <class T1, class T2>
bool operator!= (const SharedMemoryAllocator<T1>&,  const SharedMemoryAllocator<T2>&) throw() {return false;}



} // cppqedutils



#endif // CPPQEDCORE_UTILS_SHAREDMEMORYALLOCATOR_H_INCLUDED
