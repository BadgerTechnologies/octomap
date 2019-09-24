#include <octomap/OcTreeDistance.h>
#include <octomap/OcTreeDataNode.h>
#include <octomap/OcTreeKey.h>

namespace octomap {
  double
  OcTreeDistance::getMaxDistance()
  {
    return max_distance;
  }

  void
  OcTreeDistance::setMaxDistance(double max_dist)
  {
    if(max_dist > 0)
      max_distance = max_dist;
  }

  bool
  OcTreeDistance::setDistanceForKey(OcTreeKey key, double distance)
  {
    distance = std::min(std::max(distance, 0.0), max_distance);
    if(!root && !this->createRootNode())
    {
      return false;
    }
    if(!setDistanceForKeyRecursive(root, key, 0, distance))
        return false;
    return true;
  }

  bool
  OcTreeDistance::setDistanceForKeyRecursive(OcTreeDataNode<double>* node, OcTreeKey key, unsigned int depth, double distance)
  {
    if(!node)
      return false;

    if(depth < this->tree_depth)
    {
      unsigned int childPos = computeChildIdx(key, this->tree_depth-1 - depth);
      if(!this->nodeChildExists(node, childPos))
      {
        this->createNodeChild(node, childPos);
      }
      bool success = setDistanceForKeyRecursive(this->getNodeChild(node, childPos), key, depth+1, distance);
      this->pruneNode(node);
      return success;
    }
    else
    {
      node->setValue(distance);
      return true;
    }
  }
}
