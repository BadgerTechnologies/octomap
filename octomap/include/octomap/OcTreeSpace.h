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

#ifndef OCTOMAP_OCTREE_SPACE_H
#define OCTOMAP_OCTREE_SPACE_H

#include <cmath>
#include "octomap_types.h"
#include "OcTreeKey.h"

namespace octomap
{

/**
 * OcTreeSpace represents an octomap grid over space to a maximum depth.
 * Inhertied from by all OcTree varieties.
 * Having a pointer to the OcTreeSpace base of an OcTree allows code to do
 * coordinate space conversion without knowing exactly what kind of OcTree data
 * is being stored.
 */
class OcTreeSpace
{
public:
  OcTreeSpace(double res, unsigned int depth = KEY_BIT_WIDTH);
  virtual ~OcTreeSpace()
  {
  }

  OcTreeSpace(const OcTreeSpace& rhs);

  inline bool operator==(const OcTreeSpace& other) const
  {
    return tree_depth == other.tree_depth && resolution == other.resolution;
  }

  /// Change the resolution of the OcTreeSpace, scaling all voxels.
  /// This will not preserve the (metric) scale!
  virtual void setResolution(double r);
  virtual double getResolution() const
  {
    return resolution;
  }

  virtual unsigned int getTreeDepth() const
  {
    return tree_depth;
  }
  virtual void setTreeDepth(unsigned int depth);

  virtual key_type getCenterKey() const { return tree_max_val; }

  inline double getNodeSize(unsigned depth) const
  {
    if (depth > tree_depth)
      depth = tree_depth;
    return sizeLookupTable[depth];
  }

  // Always returns the key_type that would go with coordinate.
  // The return value will be aligned to the correct depth if depth is non-negative.
  // The return value will be wrapped if in_bounds is false.
  // The other coordToKey methods build on this basic function
  inline key_type coordToKeyAtDepthBoundsCheck(double coordinate, int depth = -1, bool* in_bounds = nullptr) const
  {
    assert(depth == -1 || depth <= (int)tree_depth);
    double tree_center = static_cast<double>(tree_max_val);
    double scaled_coord = std::floor(resolution_factor * coordinate) + tree_center;
    if (in_bounds)
      *in_bounds = (scaled_coord >= 0.0 && scaled_coord < 2.0 * tree_center);
    key_type keyval = static_cast<key_type>(scaled_coord);
    if (depth >= 0 && depth < (int)tree_depth)
    {
      keyval = adjustKeyAtDepth(keyval, depth);
    }
    return keyval;
  }

  /// Converts from a single coordinate into a discrete key
  inline key_type coordToKey(double coordinate) const
  {
    return coordToKeyAtDepthBoundsCheck(coordinate);
  }

  /// Converts from a single coordinate into a discrete key at a given depth
  inline key_type coordToKey(double coordinate, unsigned depth) const
  {
    return coordToKeyAtDepthBoundsCheck(coordinate, depth);
  }


  /// Converts from a 3D coordinate into a 3D addressing key
  inline OcTreeKey coordToKey(const point3d& coord) const
  {
    return OcTreeKey(coordToKey(coord(0)), coordToKey(coord(1)), coordToKey(coord(2)));
  }

  /// Converts from a 3D coordinate into a 3D addressing key
  inline OcTreeKey coordToKey(double x, double y, double z) const
  {
    return OcTreeKey(coordToKey(x), coordToKey(y), coordToKey(z));
  }

  /// Converts from a 3D coordinate into a 3D addressing key at a given depth
  inline OcTreeKey coordToKey(const point3d& coord, unsigned depth) const
  {
    if (depth == tree_depth)
      return coordToKey(coord);
    else
      return OcTreeKey(coordToKey(coord(0), depth), coordToKey(coord(1), depth), coordToKey(coord(2), depth));
  }

  /// Converts from a 3D coordinate into a 3D addressing key at a given depth
  inline OcTreeKey coordToKey(double x, double y, double z, unsigned depth) const
  {
    if (depth == tree_depth)
      return coordToKey(x, y, z);
    else
      return OcTreeKey(coordToKey(x, depth), coordToKey(y, depth), coordToKey(z, depth));
  }

  /**
   * Adjusts a 3D key from the lowest level to correspond to a higher depth (by
   * shifting the key values)
   *
   * @param key Input key, at the lowest tree level
   * @param depth Target depth level for the new key
   * @return Key for the new depth level
   */
  inline OcTreeKey adjustKeyAtDepth(const OcTreeKey& key, unsigned int depth) const
  {
    if (depth >= tree_depth)
      return key;

    return OcTreeKey(adjustKeyAtDepth(key[0], depth), adjustKeyAtDepth(key[1], depth), adjustKeyAtDepth(key[2], depth));
  }

  /**
   * Adjusts a single key value from the lowest level to correspond to a higher depth (by
   * shifting the key value)
   *
   * @param key Input key, at the lowest tree level
   * @param depth Target depth level for the new key
   * @return Key for the new depth level
   */
  inline key_type adjustKeyAtDepth(key_type key, unsigned int depth) const
  {
    unsigned int diff = tree_depth - depth;

    if (diff == 0)
      return key;
    if (diff >= tree_depth)
      return tree_max_val;
    else
      return (((key - tree_max_val) >> diff) << diff) + (1 << (diff - 1)) + tree_max_val;
  }

  /**
   * Converts a 3D coordinate into a 3D OcTreeKey, with boundary checking.
   *
   * @param point 3d coordinate of a point
   * @param key values that will be computed, an array of fixed size 3.
   * @return true if point is within the octree (valid), false otherwise
   */
  inline bool coordToKeyChecked(const point3d& point, OcTreeKey& key) const
  {
    for (unsigned int i = 0; i < 3; i++)
    {
      if (!coordToKeyChecked(point(i), key[i]))
        return false;
    }
    return true;
  }

  /**
   * Converts a 3D coordinate into a 3D OcTreeKey at a certain depth, with boundary checking.
   *
   * @param point 3d coordinate of a point
   * @param depth level of the key from the top
   * @param key values that will be computed, an array of fixed size 3.
   * @return true if point is within the octree (valid), false otherwise
   */
  inline bool coordToKeyChecked(const point3d& point, unsigned depth, OcTreeKey& key) const
  {
    for (unsigned int i = 0; i < 3; i++)
    {
      if (!coordToKeyChecked(point(i), depth, key[i]))
        return false;
    }
    return true;
  }

  /**
   * Converts a 3D coordinate into a 3D OcTreeKey, with boundary checking.
   *
   * @param x
   * @param y
   * @param z
   * @param key values that will be computed, an array of fixed size 3.
   * @return true if point is within the octree (valid), false otherwise
   */
  inline bool coordToKeyChecked(double x, double y, double z, OcTreeKey& key) const
  {
    if (!(coordToKeyChecked(x, key[0]) && coordToKeyChecked(y, key[1]) && coordToKeyChecked(z, key[2])))
    {
      return false;
    }
    else
    {
      return true;
    }
  }

  /**
   * Converts a 3D coordinate into a 3D OcTreeKey at a certain depth, with boundary checking.
   *
   * @param x
   * @param y
   * @param z
   * @param depth level of the key from the top
   * @param key values that will be computed, an array of fixed size 3.
   * @return true if point is within the octree (valid), false otherwise
   */
  inline bool coordToKeyChecked(double x, double y, double z, unsigned depth, OcTreeKey& key) const
  {
    if (!(coordToKeyChecked(x, depth, key[0]) && coordToKeyChecked(y, depth, key[1]) &&
          coordToKeyChecked(z, depth, key[2])))
    {
      return false;
    }
    else
    {
      return true;
    }
  }

  /**
   * Converts a single coordinate into a discrete addressing key, with boundary checking.
   *
   * @param coordinate 3d coordinate of a point
   * @param key discrete adressing key, result
   * @return true if coordinate is within the octree bounds (valid), false otherwise
   */
  inline bool coordToKeyChecked(double coordinate, key_type& keyval) const
  {
    bool rv;
    key_type new_keyval = coordToKeyAtDepthBoundsCheck(coordinate, -1, &rv);
    if (rv)
    {
      keyval = new_keyval;
    }
    return rv;
  }

  /**
   * Converts a single coordinate into a discrete addressing key, with boundary checking.
   *
   * @param coordinate 3d coordinate of a point
   * @param depth level of the key from the top
   * @param key discrete adressing key, result
   * @return true if coordinate is within the octree bounds (valid), false otherwise
   */
  inline bool coordToKeyChecked(double coordinate, unsigned depth, key_type& keyval) const
  {
    bool rv;
    key_type new_keyval = coordToKeyAtDepthBoundsCheck(coordinate, depth, &rv);
    if (rv)
    {
      keyval = new_keyval;
    }
    return rv;
  }

  /**
   * Converts a 3D coordinate into a 3D OcTreeKey, clamping to the boundary.
   *
   * @param point 3d coordinate of a point
   * @param key values that will be computed, an array of fixed size 3.
   * @return true if point is within the octree (valid), false otherwise
   */
  inline void coordToKeyClamped(const point3d& point, OcTreeKey& key) const
  {
    for (unsigned int i = 0; i < 3; i++)
    {
      coordToKeyClamped(point(i), key[i]);
    }
  }

  /**
   * Converts a 3D coordinate into a 3D OcTreeKey at a certain depth, clamping to the boundary.
   *
   * @param point 3d coordinate of a point
   * @param depth level of the key from the top
   * @param key values that will be computed, an array of fixed size 3.
   * @return true if point is within the octree (valid), false otherwise
   */
  inline void coordToKeyClamped(const point3d& point, unsigned depth, OcTreeKey& key) const
  {
    for (unsigned int i = 0; i < 3; i++)
    {
      coordToKeyClamped(point(i), depth, key[i]);
    }
  }

  /**
   * Converts a 3D coordinate into a 3D OcTreeKey, clamping to the boundary.
   *
   * @param x
   * @param y
   * @param z
   * @param key values that will be computed, an array of fixed size 3.
   * @return true if point is within the octree (valid), false otherwise
   */
  inline void coordToKeyClamped(double x, double y, double z, OcTreeKey& key) const
  {
    coordToKeyClamped(x, key[0]);
    coordToKeyClamped(y, key[1]);
    coordToKeyClamped(z, key[2]);
  }

  /**
   * Converts a 3D coordinate into a 3D OcTreeKey at a certain depth, clamping to the boundary.
   *
   * @param x
   * @param y
   * @param z
   * @param depth level of the key from the top
   * @param key values that will be computed, an array of fixed size 3.
   * @return true if point is within the octree (valid), false otherwise
   */
  inline void coordToKeyClamped(double x, double y, double z, unsigned depth, OcTreeKey& key) const
  {
    coordToKeyClamped(x, depth, key[0]);
    coordToKeyClamped(y, depth, key[1]);
    coordToKeyClamped(z, depth, key[2]);
  }

  /**
   * Converts a single coordinate into a discrete addressing key, clamping to the boundary.
   *
   * @param coordinate 3d coordinate of a point
   * @param key discrete adressing key, result
   * @return true if coordinate is within the octree bounds (valid), false otherwise
   */
  inline void coordToKeyClamped(double coordinate, key_type& keyval) const
  {
    double min = keyToCoord(0);
    double max = keyToCoord((tree_max_val - 1) + tree_max_val);

    if (coordinate > max)
    {
      coordinate = max;
    }
    else if (coordinate < min)
    {
      coordinate = min;
    }

    keyval = coordToKey(coordinate);
  }

  /**
   * Converts a single coordinate into a discrete addressing key, clamping to the boundary.
   *
   * @param coordinate 3d coordinate of a point
   * @param depth level of the key from the top
   * @param key discrete adressing key, result
   * @return true if coordinate is within the octree bounds (valid), false otherwise
   */
  inline void coordToKeyClamped(double coordinate, unsigned depth, key_type& keyval) const
  {
    coordToKeyClamped(coordinate, keyval);
    if (depth < tree_depth)
    {
      keyval = adjustKeyAtDepth(keyval, depth);
    }
  }

  /// converts from a discrete key at a given depth into a coordinate
  /// corresponding to the key's center
  inline double keyToCoord(key_type key, unsigned depth) const
  {
    assert(depth <= tree_depth);

    // root is centered on 0 = 0.0
    if (depth == 0)
    {
      return 0.0;
    }
    else if (depth == tree_depth)
    {
      return keyToCoord(key);
    }
    else
    {
      return (floor((double(key) - double(this->tree_max_val))/ double(1 << (tree_depth - depth)))
              + 0.5) * this->getNodeSize(depth);
    }
  }

  /// converts from a discrete key at the lowest tree level into a coordinate
  /// corresponding to the key's center
  inline double keyToCoord(key_type key) const
  {
    return (double(key) - double(this->tree_max_val) + 0.5) * this->resolution;
  }

  /// converts from an addressing key at the lowest tree level into a coordinate
  /// corresponding to the key's center
  inline point3d keyToCoord(const OcTreeKey& key) const
  {
    return point3d(keyToCoord(key[0]), keyToCoord(key[1]), keyToCoord(key[2]));
  }

  /// converts from an addressing key at a given depth into a coordinate
  /// corresponding to the key's center
  inline point3d keyToCoord(const OcTreeKey& key, unsigned depth) const
  {
    return point3d(keyToCoord(key[0], depth), keyToCoord(key[1], depth), keyToCoord(key[2], depth));
  }

protected:
  unsigned int tree_depth;
  key_type tree_max_val;     ///< really center value, derived from tree_depth, for convenience
  double resolution;         ///< in meters
  double resolution_factor;  ///< = 1. / resolution

  /// contains the size of a voxel at tree level i (0: root node). tree_depth+1 levels (incl. 0)
  std::vector<double> sizeLookupTable;
};

}  // namespace octomap

#endif  // OCTOMAP_OCTREE_SPACE_H
