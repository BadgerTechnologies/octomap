/*
 * OctoMap - An Efficient Probabilistic 3D Mapping Framework Based on Octrees
 * http://octomap.github.com/
 *
 * Copyright (c) 2009-2013, K.M. Wurm and A. Hornung, University of Freiburg
 * All rights reserved.
 * License: New BSD
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * Neither the name of the University of Freiburg nor the names of its
 *       contributors may be used to endorse or promote products derived from
 *       this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */

#ifndef OCTOMAP_OCTREE_KEY_H
#define OCTOMAP_OCTREE_KEY_H

/* According to c++ standard including this header has no practical effect
 * but it can be used to determine the c++ standard library implementation.
 */
#include <ciso646>

#include <assert.h>

/* Libc++ does not implement the TR1 namespace, all c++11 related functionality
 * is instead implemented in the std namespace.
 */
#if defined(__GNUC__) && ! defined(_LIBCPP_VERSION)
  #include <tr1/unordered_set>
  #include <tr1/unordered_map>
  namespace octomap {
    namespace unordered_ns = std::tr1;
  }
#else
  #include <unordered_set>
  #include <unordered_map>
  namespace octomap {
    namespace unordered_ns = std;
  }
#endif

namespace octomap {

  using key_type = uint32_t;
  static constexpr unsigned int KEY_BIT_WIDTH = sizeof(key_type)*8;
  // Center of the key space, used as the origin of tree coordinates
  static constexpr key_type KEY_CENTER = (1 << (KEY_BIT_WIDTH-1));

  
  /**
   * OcTreeKey is a container class for internal key addressing. The keys count the
   * number of cells (voxels) from the origin as discrete address of a voxel.
   * @see OcTreeBaseImpl::coordToKey() and OcTreeBaseImpl::keyToCoord() for conversions.
   */
  class OcTreeKey {
    
  public:  
    OcTreeKey () {
      // Avoid uninitialized warning for k[3] by setting all values to zero
      for (unsigned int i=0; i<4; ++i) {
        k[i] = 0;
      }
    }
    OcTreeKey (key_type a, key_type b, key_type c){ 
      k[0] = a; 
      k[1] = b; 
      k[2] = c;         
    }
    
    OcTreeKey(const OcTreeKey& other){
      for (unsigned int i=0; i<4; ++i) {
        k[i] = other.k[i];
      }
    }
    
    bool operator== (const OcTreeKey &other) const { 
      return ((k[0] == other.k[0]) && (k[1] == other.k[1]) && (k[2] == other.k[2]));
    }
    
    bool operator!= (const OcTreeKey& other) const {
      return( (k[0] != other.k[0]) || (k[1] != other.k[1]) || (k[2] != other.k[2]));
    }
    
    OcTreeKey& operator=(const OcTreeKey& other){
      for (unsigned int i=0; i<4; ++i) {
        k[i] = other.k[i];
      }
      return *this;
    }
    
    const key_type& operator[] (unsigned int i) const { 
      return k[i];
    }
    
    key_type& operator[] (unsigned int i) { 
      return k[i];
    }

    alignas(16) key_type k[4];

    /// Provides a hash function on Keys
    struct KeyHash{
      size_t operator()(const OcTreeKey& key) const{
        // a simple hashing function 
	// explicit casts to size_t to operate on the complete range
	// constanst will be promoted according to C++ standard
        return static_cast<size_t>(key.k[0])
          + 1447*static_cast<size_t>(key.k[1])
          + 345637*static_cast<size_t>(key.k[2]);
      }
    };
    
  };
  
  /**
   * Data structure to efficiently compute the nodes to update from a scan
   * insertion using a hash set.
   * @note you need to use boost::unordered_set instead if your compiler does not
   * yet support tr1!
   */
  typedef unordered_ns::unordered_set<OcTreeKey, OcTreeKey::KeyHash> KeySet;

  /**
   * Data structrure to efficiently track changed nodes as a combination of
   * OcTreeKeys and a bool flag (to denote newly created nodes)
   *
   */
  typedef unordered_ns::unordered_map<OcTreeKey, bool, OcTreeKey::KeyHash> KeyBoolMap;


  class KeyRay {
  public:
    
    KeyRay () {
      ray.resize(maxSize);
      reset();
    }
    
    KeyRay(const KeyRay& other){
      ray = other.ray;
      size_t dSize = other.end() - other.begin();
      end_of_ray = ray.begin() + dSize;
    }
    
    void reset() {
      end_of_ray = begin();
    }
    
    void addKey(const OcTreeKey& k) {
      assert(end_of_ray != ray.end());
      *end_of_ray = k;
      ++end_of_ray;
    }

    size_t size() const { return end_of_ray - ray.begin(); }
    size_t sizeMax() const { return maxSize; }

    typedef std::vector<OcTreeKey>::iterator iterator;
    typedef std::vector<OcTreeKey>::const_iterator const_iterator;
    typedef std::vector<OcTreeKey>::reverse_iterator reverse_iterator;
    
    iterator begin() { return ray.begin(); }
    iterator end() { return end_of_ray; }
    const_iterator begin() const { return ray.begin(); }
    const_iterator end() const   { return end_of_ray; }

    reverse_iterator rbegin() { return (reverse_iterator) end_of_ray; }
    reverse_iterator rend() { return ray.rend(); }

  private:
    std::vector<OcTreeKey> ray;
    std::vector<OcTreeKey>::iterator end_of_ray;
    const static size_t maxSize = 100000;
  };

  /**
   * Computes the center offset key at given depth
   *
   * @param[in] depth of center offset key to find
   * @param[in] tree_center_key center key of the tree (also erroneously called tree_max_val)
   */
  inline key_type computeCenterOffsetKey(unsigned int depth, key_type center_key) {
    if (depth + 1 >= KEY_BIT_WIDTH) {
      // It is undefined in C++ to shift by any amount >= to the bit width of
      // the type. Return zero in this case.
      return 0;
    }
    return center_key >> (depth + 1);
  }

  /**
   * Computes the key of a child node while traversing the octree, given
   * child index and current key
   *
   * @param[in] pos index of child node (0..7)
   * @param[in] center_offset_key constant offset of octree keys
   * @param[in] parent_key current (parent) key
   * @param[out] child_key  computed child key
   */
  inline void computeChildKey (unsigned int pos, key_type center_offset_key,
                                          const OcTreeKey& parent_key, OcTreeKey& child_key) {
    const key_type offset_keys[2] = {0 - center_offset_key - (center_offset_key ? 0 : 1), center_offset_key};
    // x-axis
    child_key[0] = parent_key[0] + offset_keys[(pos & 1) != 0];
    // y-axis
    child_key[1] = parent_key[1] + offset_keys[(pos & 2) != 0];
    // z-axis
    child_key[2] = parent_key[2] + offset_keys[(pos & 4) != 0];
  }
  
  /// generate child index (between 0 and 7) from key at given tree depth
  inline uint8_t computeChildIdx(const OcTreeKey& key, int depth){
    const key_type bit = (1 << depth);
    return ((key.k[0] & bit) != 0)
         | (((key.k[1] & bit) != 0)<<1)
         | (((key.k[2] & bit) != 0)<<2);
  }

  /**
   * Generates a unique key for all keys on a certain level of the tree
   *
   * @param level from the bottom (= tree_depth - depth of key)
   * @param key input indexing key (at lowest resolution / level)
   * @return key corresponding to the input key at the given level
   */
  inline OcTreeKey computeIndexKey(key_type level, const OcTreeKey& key) {
    if (level == 0)
      return key;
    else if (level >= KEY_BIT_WIDTH)
      // avoid undefined behavior with the bit shift below
      return OcTreeKey(0,0,0);
    else {
      const key_type mask = std::numeric_limits<key_type>::max() << level;
      OcTreeKey result;
      result[0] = key[0] & mask;
      result[1] = key[1] & mask;
      result[2] = key[2] & mask;
      return result;
    }
  }

} // namespace

#endif
