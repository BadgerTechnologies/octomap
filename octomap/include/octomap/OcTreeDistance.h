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

#ifndef OCTOMAP_OCTREE_DISTANCE_H
#define OCTOMAP_OCTREE_DISTANCE_H


#include "OcTreeBaseImpl.h"
#include "AbstractOcTree.h"
#include "OcTreeDataNode.h"
#include "OcTreeKey.h"


namespace octomap {
  class OcTreeDistance : public OcTreeBaseImpl<OcTreeDataNode<double>,AbstractOcTree> {
  public:
    OcTreeDistance(double res, double max_dist) : OcTreeBaseImpl<OcTreeDataNode<double>,AbstractOcTree>(res)
    {
      max_distance = max_dist;
    };

    /// virtual constructor: creates a new object of same type
    /// (Covariant return type requires an up-to-date compiler)
    OcTreeDistance* create() const {return new OcTreeDistance(this->resolution, this->max_distance); }
    std::string getTreeType() const {return "OcTreeDistance";}
    
    double getMaxDistance();
    void setMaxDistance(double max_dist);

    bool setDistanceForKey(OcTreeKey key, double distance);
  private:
    double max_distance;
    bool setDistanceForKeyRecursive(OcTreeDataNode<double>* node, OcTreeKey key, unsigned int depth, double distance);
  };

  };


#endif
