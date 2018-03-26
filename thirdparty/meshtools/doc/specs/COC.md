CoC concept
-----------

A *Container of Container* (CoC) must define the following interface:
* `value_type` is the type of stored containers it must itself define at least:
  * `value_type` is the type of stored elements
  * `begin()` and `end()` methods that return iterators that point respectively to the first and past element of the container
  * `size() const` shall return the number of elements
* `begin()` and `end()` methods shall return iterators that point respectively to the first and past last container
* `operator[](const int i)` shall return the i<sup>th</sup> container 
* `size() const` shall return the number of containers 

FSCoC
-----
A *Fixed Sized CoC* (FSCoC) is just a vector of fixed size array. Its implementation is
straightforward using STL classes:
```cpp
#include <vector>
#include <array>
template <typename T, std::size_t N>
using FSCoC = std::vector< std::array<T,N> >;
```
It has the CoC interface by construction.

