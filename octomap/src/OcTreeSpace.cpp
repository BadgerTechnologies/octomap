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

#include <octomap/OcTreeSpace.h>

namespace octomap
{

OcTreeSpace::OcTreeSpace(double res, unsigned int depth)
{
  tree_depth = depth;
  tree_max_val = (1U << (depth - 1));
  setResolution(res);
}

OcTreeSpace::OcTreeSpace(const OcTreeSpace& rhs)
  : tree_depth(rhs.tree_depth)
  , tree_max_val(rhs.tree_max_val)
  , resolution(rhs.resolution)
  , resolution_factor(rhs.resolution_factor)
  , sizeLookupTable(rhs.sizeLookupTable)
{
}

void OcTreeSpace::setResolution(double r)
{
  resolution = r;
  resolution_factor = 1. / resolution;

  // init node size lookup table:
  sizeLookupTable.resize(tree_depth + 1);
  for (unsigned i = 1; i <= tree_depth; ++i)
  {
    sizeLookupTable[i] = resolution * double(((size_t)1) << (tree_depth - i));
  }
  // to keep from getting the wrong answer when at depth 0, double the
  // previous result instead of possibly shifting the 1 off the end
  sizeLookupTable[0] = 2.0 * sizeLookupTable[1];
}

void OcTreeSpace::setTreeDepth(unsigned int depth)
{
  if (tree_depth == depth)
  {
    // Nothing to do
    return;
  }

  tree_depth = depth;
  tree_max_val = (1U << (depth - 1));
  // need to recalculate lut
  setResolution(resolution);
}

}  // namespace octomap
