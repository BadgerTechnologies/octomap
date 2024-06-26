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

#include <bitset>
#include <algorithm>

#include <octomap/MCTables.h>

namespace octomap {

  template <class NODE>
  OccupancyOcTreeBase<NODE>::OccupancyOcTreeBase(double in_resolution)
    : OcTreeBaseImpl<NODE,AbstractOccupancyOcTree>(in_resolution), use_bbx_limit(false), use_change_detection(false)
  {

  }

  template <class NODE>
  OccupancyOcTreeBase<NODE>::OccupancyOcTreeBase(double in_resolution, unsigned int in_tree_depth, unsigned int in_tree_max_val)
    : OcTreeBaseImpl<NODE,AbstractOccupancyOcTree>(in_resolution, in_tree_depth, in_tree_max_val), use_bbx_limit(false), use_change_detection(false)
  {

  }

  template <class NODE>
  OccupancyOcTreeBase<NODE>::~OccupancyOcTreeBase(){
  }

  template <class NODE>
  OccupancyOcTreeBase<NODE>::OccupancyOcTreeBase(const OccupancyOcTreeBase<NODE>& rhs) :
  OcTreeBaseImpl<NODE,AbstractOccupancyOcTree>(rhs), use_bbx_limit(rhs.use_bbx_limit),
    bbx_min(rhs.bbx_min), bbx_max(rhs.bbx_max),
    bbx_min_key(rhs.bbx_min_key), bbx_max_key(rhs.bbx_max_key),
    use_change_detection(rhs.use_change_detection), changed_keys(rhs.changed_keys)
  {
    this->clamping_thres_min = rhs.clamping_thres_min;
    this->clamping_thres_max = rhs.clamping_thres_max;
    this->prob_hit_log = rhs.prob_hit_log;
    this->prob_miss_log = rhs.prob_miss_log;
    this->occ_prob_thres_log = rhs.occ_prob_thres_log;


  }

  template <class NODE>
  void OccupancyOcTreeBase<NODE>::insertPointCloud(const ScanNode& scan, double maxrange, bool lazy_eval, bool discretize) {
    // performs transformation to data and sensor origin first
    Pointcloud& cloud = *(scan.scan);
    pose6d frame_origin = scan.pose;
    point3d sensor_origin = frame_origin.inv().transform(scan.pose.trans());
    insertPointCloud(cloud, sensor_origin, frame_origin, maxrange, lazy_eval, discretize);
  }


  template <class NODE>
  void OccupancyOcTreeBase<NODE>::insertPointCloud(const Pointcloud& scan, const octomap::point3d& sensor_origin,
                                             double maxrange, bool lazy_eval, bool discretize) {

    KeySet free_cells, occupied_cells;
    if (discretize)
      computeDiscreteUpdate(scan, sensor_origin, free_cells, occupied_cells, maxrange);
    else
      computeUpdate(scan, sensor_origin, free_cells, occupied_cells, maxrange);

    // insert data into tree  -----------------------
    for (KeySet::iterator it = free_cells.begin(); it != free_cells.end(); ++it) {
      updateNode(*it, false, lazy_eval);
    }
    for (KeySet::iterator it = occupied_cells.begin(); it != occupied_cells.end(); ++it) {
      updateNode(*it, true, lazy_eval);
    }
  }

  template <class NODE>
  void OccupancyOcTreeBase<NODE>::insertPointCloud(const Pointcloud& pc, const point3d& sensor_origin, const pose6d& frame_origin,
                                             double maxrange, bool lazy_eval, bool discretize) {
    // performs transformation to data and sensor origin first
    Pointcloud transformed_scan (pc);
    transformed_scan.transform(frame_origin);
    point3d transformed_sensor_origin = frame_origin.transform(sensor_origin);
    insertPointCloud(transformed_scan, transformed_sensor_origin, maxrange, lazy_eval, discretize);
  }


  template <class NODE>
  void OccupancyOcTreeBase<NODE>::insertPointCloudRays(const Pointcloud& pc, const point3d& origin, double /* maxrange */, bool lazy_eval) {
    if (pc.size() < 1)
      return;

#ifdef _OPENMP
    omp_set_num_threads(this->keyrays.size());
    #pragma omp parallel for
#endif
    for (int i = 0; i < (int)pc.size(); ++i) {
      const point3d& p = pc[i];
      unsigned threadIdx = 0;
#ifdef _OPENMP
      threadIdx = omp_get_thread_num();
#endif
      KeyRay* keyray = &(this->keyrays.at(threadIdx));

      if (this->computeRayKeys(origin, p, *keyray)){
#ifdef _OPENMP
        #pragma omp critical
#endif
        {
          for(KeyRay::iterator it=keyray->begin(); it != keyray->end(); it++) {
            updateNode(*it, false, lazy_eval); // insert freespace measurement
          }
          updateNode(p, true, lazy_eval); // update endpoint to be occupied
        }
      }

    }
  }

  template <class NODE>
  void OccupancyOcTreeBase<NODE>::computeDiscreteUpdate(const Pointcloud& scan, const octomap::point3d& origin,
                                                KeySet& free_cells, KeySet& occupied_cells,
                                                double maxrange)
 {
   Pointcloud discretePC;
   discretePC.reserve(scan.size());
   KeySet endpoints;

   for (int i = 0; i < (int)scan.size(); ++i) {
     OcTreeKey k = this->coordToKey(scan[i]);
     std::pair<KeySet::iterator,bool> ret = endpoints.insert(k);
     if (ret.second){ // insertion took place => k was not in set
       discretePC.push_back(this->keyToCoord(k));
     }
   }

   computeUpdate(discretePC, origin, free_cells, occupied_cells, maxrange);
 }


  template <class NODE>
  void OccupancyOcTreeBase<NODE>::computeUpdate(const Pointcloud& scan, const octomap::point3d& origin,
                                                KeySet& free_cells, KeySet& occupied_cells,
                                                double maxrange)
  {



#ifdef _OPENMP
    omp_set_num_threads(this->keyrays.size());
    #pragma omp parallel for schedule(guided)
#endif
    for (int i = 0; i < (int)scan.size(); ++i) {
      const point3d& p = scan[i];
      unsigned threadIdx = 0;
#ifdef _OPENMP
      threadIdx = omp_get_thread_num();
#endif
      KeyRay* keyray = &(this->keyrays.at(threadIdx));


      if (!use_bbx_limit) { // no BBX specified
        if ((maxrange < 0.0) || ((p - origin).norm() <= maxrange) ) { // is not maxrange meas.
          // free cells
          if (this->computeRayKeys(origin, p, *keyray)){
#ifdef _OPENMP
            #pragma omp critical (free_insert)
#endif
            {
              free_cells.insert(keyray->begin(), keyray->end());
            }
          }
          // occupied endpoint
          OcTreeKey key;
          if (this->coordToKeyChecked(p, key)){
#ifdef _OPENMP
            #pragma omp critical (occupied_insert)
#endif
            {
              occupied_cells.insert(key);
            }
          }
        } else { // user set a maxrange and length is above
          point3d direction = (p - origin).normalized ();
          point3d new_end = origin + direction * maxrange;
          if (this->computeRayKeys(origin, new_end, *keyray)){
#ifdef _OPENMP
            #pragma omp critical (free_insert)
#endif
            {
              free_cells.insert(keyray->begin(), keyray->end());
            }
          }
        } // end if maxrange
      } else { // BBX was set
        // endpoint in bbx and not maxrange?
        if ( inBBX(p) && ((maxrange < 0.0) || ((p - origin).norm () <= maxrange) ) )  {

          // occupied endpoint
          OcTreeKey key;
          if (this->coordToKeyChecked(p, key)){
#ifdef _OPENMP
            #pragma omp critical (occupied_insert)
#endif
            {
              occupied_cells.insert(key);
            }
          }

          // update freespace, break as soon as bbx limit is reached
          if (this->computeRayKeys(origin, p, *keyray)){
            for(KeyRay::reverse_iterator rit=keyray->rbegin(); rit != keyray->rend(); rit++) {
              if (inBBX(*rit)) {
#ifdef _OPENMP
                #pragma omp critical (free_insert)
#endif
                {
                  free_cells.insert(*rit);
                }
              }
              else break;
            }
          } // end if compute ray
        } // end if in BBX and not maxrange
      } // end bbx case

    } // end for all points, end of parallel OMP loop

    // prefer occupied cells over free ones (and make sets disjunct)
    for(KeySet::iterator it = free_cells.begin(), end=free_cells.end(); it!= end; ){
      if (occupied_cells.find(*it) != occupied_cells.end()){
        it = free_cells.erase(it);
      } else {
        ++it;
      }
    }
  }

  template <class NODE>
  NODE* OccupancyOcTreeBase<NODE>::setNodeValue(const OcTreeKey& key, float log_odds_value, bool lazy_eval) {
    // clamp log odds within range:
    log_odds_value = std::min(std::max(log_odds_value, this->clamping_thres_min), this->clamping_thres_max);

    bool createdRoot = false;
    if (this->root == NULL){
      this->root = this->allocNode();
      this->tree_size++;
      createdRoot = true;
    }

    return setNodeValueRecurs(this->root, createdRoot, key, 0, log_odds_value, lazy_eval);
  }

  template <class NODE>
  NODE* OccupancyOcTreeBase<NODE>::setNodeValue(const point3d& value, float log_odds_value, bool lazy_eval) {
    OcTreeKey key;
    if (!this->coordToKeyChecked(value, key))
      return NULL;

    return setNodeValue(key, log_odds_value, lazy_eval);
  }

  template <class NODE>
  NODE* OccupancyOcTreeBase<NODE>::setNodeValue(double x, double y, double z, float log_odds_value, bool lazy_eval) {
    OcTreeKey key;
    if (!this->coordToKeyChecked(x, y, z, key))
      return NULL;

    return setNodeValue(key, log_odds_value, lazy_eval);
  }

  template <class NODE>
  NODE* OccupancyOcTreeBase<NODE>::setNodeValueAtDepth(const OcTreeKey& key, unsigned int depth,
                                                       float log_odds_value, bool lazy_eval) {
    // clamp log odds within range:
    log_odds_value = std::min(std::max(log_odds_value, this->clamping_thres_min), this->clamping_thres_max);

    if (depth > this->tree_depth)
      depth = this->tree_depth;

    bool createdRoot = false;
    if (this->root == NULL){
      this->root = this->allocNode();
      this->tree_size++;
      createdRoot = true;
    }

    return setNodeValueAtDepthRecurs(this->root, createdRoot, key, 0, depth, log_odds_value, lazy_eval);
  }

  template <class NODE>
  NODE* OccupancyOcTreeBase<NODE>::updateNode(const OcTreeKey& key, float log_odds_update, bool lazy_eval) {
    // early abort (no change will happen).
    // may cause an overhead in some configuration, but more often helps
    NODE* leaf = this->search(key);
    // no change: node already at threshold
    if (leaf
        && ((log_odds_update >= 0 && leaf->getLogOdds() >= this->clamping_thres_max)
        || ( log_odds_update <= 0 && leaf->getLogOdds() <= this->clamping_thres_min)))
    {
      return leaf;
    }

    bool createdRoot = false;
    if (this->root == NULL){
      this->root = this->allocNode();
      this->tree_size++;
      createdRoot = true;
    }

    return updateNodeRecurs(this->root, createdRoot, key, 0, log_odds_update, lazy_eval);
  }

  template <class NODE>
  NODE* OccupancyOcTreeBase<NODE>::updateNode(const point3d& value, float log_odds_update, bool lazy_eval) {
    OcTreeKey key;
    if (!this->coordToKeyChecked(value, key))
      return NULL;

    return updateNode(key, log_odds_update, lazy_eval);
  }

  template <class NODE>
  NODE* OccupancyOcTreeBase<NODE>::updateNode(double x, double y, double z, float log_odds_update, bool lazy_eval) {
    OcTreeKey key;
    if (!this->coordToKeyChecked(x, y, z, key))
      return NULL;

    return updateNode(key, log_odds_update, lazy_eval);
  }

  template <class NODE>
  NODE* OccupancyOcTreeBase<NODE>::updateNode(const OcTreeKey& key, bool occupied, bool lazy_eval) {
    float logOdds = this->prob_miss_log;
    if (occupied)
      logOdds = this->prob_hit_log;

    return updateNode(key, logOdds, lazy_eval);
  }

  template <class NODE>
  NODE* OccupancyOcTreeBase<NODE>::updateNode(const point3d& value, bool occupied, bool lazy_eval) {
    OcTreeKey key;
    if (!this->coordToKeyChecked(value, key))
      return NULL;
    return updateNode(key, occupied, lazy_eval);
  }

  template <class NODE>
  NODE* OccupancyOcTreeBase<NODE>::updateNode(double x, double y, double z, bool occupied, bool lazy_eval) {
    OcTreeKey key;
    if (!this->coordToKeyChecked(x, y, z, key))
      return NULL;
    return updateNode(key, occupied, lazy_eval);
  }

  template <class NODE>
  NODE* OccupancyOcTreeBase<NODE>::updateNodeRecurs(NODE* node, bool node_just_created, const OcTreeKey& key,
                                                    unsigned int depth, const float& log_odds_update, bool lazy_eval) {
    bool created_node = false;

    assert(node);

    // follow down to last level
    if (depth < this->tree_depth) {
      unsigned int pos = computeChildIdx(key, this->tree_depth -1 - depth);
      if (!this->nodeChildExists(node, pos)) {
        // child does not exist, but maybe it's a pruned node?
        if (!this->nodeHasChildren(node) && !node_just_created ) {
          // current node does not have children AND it is not a new node
          // -> expand pruned node
          this->expandNode(node);
        }
        else {
          // not a pruned node, create requested child
          this->createNodeChild(node, pos);
          created_node = true;
        }
      }

      if (lazy_eval)
        return updateNodeRecurs(this->getNodeChild(node, pos), created_node, key, depth+1, log_odds_update, lazy_eval);
      else {
        NODE* retval = updateNodeRecurs(this->getNodeChild(node, pos), created_node, key, depth+1, log_odds_update, lazy_eval);
        // prune node if possible, otherwise set own probability
        // note: combining both did not lead to a speedup!
        if (this->pruneNode(node)){
          // return pointer to current parent (pruned), the just updated node no longer exists
          retval = node;
        } else{
          node->updateOccupancyChildren();
        }

        return retval;
      }
    }

    // at last level, update node, end of recursion
    else {
      updateNodeLogOddsAndTrackChanges(node, log_odds_update, node_just_created, key, depth);
      return node;
    }
  }

  // TODO: mostly copy of updateNodeRecurs => merge code or general tree modifier / traversal
  template <class NODE>
  NODE* OccupancyOcTreeBase<NODE>::setNodeValueRecurs(NODE* node, bool node_just_created, const OcTreeKey& key,
                                                    unsigned int depth, const float& log_odds_value, bool lazy_eval) {
    bool created_node = false;

    assert(node);

    // follow down to last level
    if (depth < this->tree_depth) {
      unsigned int pos = computeChildIdx(key, this->tree_depth -1 - depth);
      if (!this->nodeChildExists(node, pos)) {
        // child does not exist, but maybe it's a pruned node?
        if (!this->nodeHasChildren(node) && !node_just_created ) {
          // current node does not have children AND it is not a new node
          // -> expand pruned node
          this->expandNode(node);
        }
        else {
          // not a pruned node, create requested child
          this->createNodeChild(node, pos);
          created_node = true;
        }
      }

      if (lazy_eval)
        return setNodeValueRecurs(this->getNodeChild(node, pos), created_node, key, depth+1, log_odds_value, lazy_eval);
      else {
        NODE* retval = setNodeValueRecurs(this->getNodeChild(node, pos), created_node, key, depth+1, log_odds_value, lazy_eval);
        // prune node if possible, otherwise set own probability
        // note: combining both did not lead to a speedup!
        if (this->pruneNode(node)){
          // return pointer to current parent (pruned), the just updated node no longer exists
          retval = node;
        } else{
          node->updateOccupancyChildren();
        }

        return retval;
      }
    }

    // at last level, update node, end of recursion
    else {
      if (use_change_detection || nodeValueChangeCallback) {
        bool occBefore = this->isNodeOccupied(node);
        float valBefore = node->getLogOdds();
        node->setLogOdds(log_odds_value);

        if (node->getLogOdds() != valBefore || node_just_created)
        {
          valueChangeCallbackWrapper(key, depth, node_just_created, valBefore, occBefore, node->getLogOdds(), this->isNodeOccupied(node));
        }

        if (use_change_detection) {
          if (node_just_created){  // new node
            changed_keys.insert(std::pair<OcTreeKey,bool>(key, true));
          } else if (occBefore != this->isNodeOccupied(node)) {  // occupancy changed, track it
            KeyBoolMap::iterator it = changed_keys.find(key);
            if (it == changed_keys.end())
              changed_keys.insert(std::pair<OcTreeKey,bool>(key, false));
            else if (it->second == false)
              changed_keys.erase(it);
          }
        }
      } else {
        node->setLogOdds(log_odds_value);
      }
      return node;
    }
  }

  template <class NODE>
  NODE* OccupancyOcTreeBase<NODE>::setNodeValueAtDepthRecurs(NODE* node,
                                                             bool node_just_created,
                                                             const OcTreeKey& key,
                                                             unsigned int current_depth,
                                                             unsigned int target_depth,
                                                             float log_odds_value,
                                                             bool lazy_eval) {
    bool created_node = false;
    if (current_depth == target_depth || current_depth == this->tree_depth) {
      // At the target node (or end of tree). Delete all children and set the log odds
      this->deleteNodeChildren(node);
      node->setLogOdds(log_odds_value);
      return node;
    }

    // Not yet at the target node. Decend the tree.
    unsigned int pos = computeChildIdx(key, this->tree_depth - 1 - current_depth);
    if (!this->nodeChildExists(node, pos)) {
      if (!this->nodeHasChildren(node) && !node_just_created) {
        // was a pruned node, expand it
        this->expandNode(node);
      } else {
        // create new child
        this->createNodeChild(node, pos);
        created_node = true;
      }
    }
    NODE* retval = setNodeValueAtDepthRecurs(this->getNodeChild(node, pos),
                                             created_node,
                                             key,
                                             current_depth + 1,
                                             target_depth,
                                             log_odds_value,
                                             lazy_eval);
    if (!lazy_eval) {
      if (this->pruneNode(node)) {
        retval = node;
      } else {
        node->updateOccupancyChildren();
      }
    }
    return retval;
  }

  template <class NODE>
  void OccupancyOcTreeBase<NODE>::setTreeValues(const OccupancyOcTreeBase<NODE>* value_tree,
                                                bool maximum_only, bool delete_first,
                                                CopyValueFunction copy_value_function){
    setTreeValues(
        static_cast<const OccupancyOcTreeBase<NODE>*>(value_tree),
        static_cast<const OccupancyOcTreeBase<NODE>*>(NULL),
        maximum_only, delete_first, copy_value_function);
  }

  template <class NODE>
  void OccupancyOcTreeBase<NODE>::setTreeValues(OccupancyOcTreeBase<NODE>* value_tree,
                                                bool maximum_only, bool delete_first,
                                                CopyValueFunction copy_value_function){
    if (value_tree && this->getTreeDepth() != value_tree->getTreeDepth()) {
      value_tree->setTreeDepth(this->getTreeDepth());
    }
    setTreeValues(
        static_cast<const OccupancyOcTreeBase<NODE>*>(value_tree),
        static_cast<const OccupancyOcTreeBase<NODE>*>(NULL),
        maximum_only, delete_first, copy_value_function);
  }


  template <class NODE>
  void OccupancyOcTreeBase<NODE>::setTreeValues(const OccupancyOcTreeBase<NODE>* value_tree,
                                                const OccupancyOcTreeBase<NODE>* bounds_tree,
                                                bool maximum_only, bool delete_first,
                                                CopyValueFunction copy_value_function){
    octomap::OcTreeKey root_key(this->tree_max_val, this->tree_max_val, this->tree_max_val);

    // setTreeValues only makes sense when the trees have nearly the same resolution
    if (value_tree && (this->getResolution() - value_tree->getResolution()) > 1e-6) {
      OCTOMAP_ERROR_STR(
          "setTreeValues: tree resolution mismatch with value tree: " <<
          this->getResolution() << " != " <<
          value_tree->getResolution());
      return;
    }
    if (bounds_tree && (this->getResolution() - bounds_tree->getResolution()) > 1e-6) {
      OCTOMAP_ERROR_STR(
          "setTreeValues: tree resolution mismatch with bounds tree: " <<
          this->getResolution() << " != " <<
          bounds_tree->getResolution());
      return;
    }

    // setTreeValues only works when the value_tree and bounds_tree have the
    // same depth as this tree.
    if (value_tree && this->getTreeDepth() != value_tree->getTreeDepth()) {
      OCTOMAP_ERROR_STR(
          "setTreeValues: tree depth mismatch with value tree: " <<
          this->getTreeDepth() << "!= " <<
          value_tree->getTreeDepth());
      return;
    }
    if (bounds_tree && this->getTreeDepth() != bounds_tree->getTreeDepth()) {
      OCTOMAP_ERROR_STR(
          "setTreeValues: tree depth mismatch with bounds tree: " <<
          this->getTreeDepth() << "!= " <<
          bounds_tree->getTreeDepth());
      return;
    }

    // delete_first implies maximum_only = false.
    if (delete_first)
      maximum_only = false;

    bool created_root = false;
    if (this->root == NULL)
    {
      this->root = this->allocNode();
      this->tree_size++;
      created_root = true;
    }

    const NODE* value_node = NULL;
    if (value_tree != NULL) {
      value_node = value_tree->root;
    }

    const NODE* bounds_node = NULL;
    if (bounds_tree != NULL) {
      bounds_node = bounds_tree->root;
    }

    this->root = setTreeValuesRecurs(this->root,
                                     created_root,
                                     root_key,
                                     0,
                                     value_tree,
                                     bounds_tree,
                                     value_node,
                                     bounds_node,
                                     maximum_only,
                                     delete_first,
                                     copy_value_function);
  }

  template <class NODE>
  void OccupancyOcTreeBase<NODE>::setTreeValues(OccupancyOcTreeBase<NODE>* value_tree,
                                                OccupancyOcTreeBase<NODE>* bounds_tree,
                                                bool maximum_only, bool delete_first,
                                                CopyValueFunction copy_value_function){
    if (value_tree && this->getTreeDepth() != value_tree->getTreeDepth()) {
      value_tree->setTreeDepth(this->getTreeDepth());
    }
    if (bounds_tree && this->getTreeDepth() != bounds_tree->getTreeDepth()) {
      bounds_tree->setTreeDepth(this->getTreeDepth());
    }
    setTreeValues(
        static_cast<const OccupancyOcTreeBase<NODE>*>(value_tree),
        static_cast<const OccupancyOcTreeBase<NODE>*>(bounds_tree),
        maximum_only, delete_first, copy_value_function);
  }

  template <class NODE>
  NODE* OccupancyOcTreeBase<NODE>::setTreeValuesRecurs(NODE* node,
                                                       bool node_just_created,
                                                       const OcTreeKey& key,
                                                       unsigned int depth,
                                                       const OccupancyOcTreeBase<NODE>* value_tree,
                                                       const OccupancyOcTreeBase<NODE>* bounds_tree,
                                                       const NODE* value_node,
                                                       const NODE* bounds_node,
                                                       bool maximum_only,
                                                       bool delete_first,
                                                       const CopyValueFunction& copy_value_function){
    if (bounds_tree != NULL && bounds_node == NULL) {
      // We are out-of-bounds.
      // Must check this first to not attempt referencing bounds_node
      // Be sure to delete a newly-created node, as nothing will be added to it.
      if (node != NULL && node_just_created) {
        this->deleteNodeRecurs(node);
        node = NULL;
      }
      return node;
    }
    if (bounds_tree != NULL && bounds_tree->nodeHasChildren(bounds_node)) {
      // We aren't yet in a leaf of the bounds tree, descend
      // First expand if our node is a pruned leaf.
      if (!this->nodeHasChildren(node) && !node_just_created) {
        // Our node is a leaf. Expand it as we will need
        // to update some of the sub-nodes.
        this->expandNode(node);
      }
      for (unsigned int i=0; i<8; ++i) {
        if (bounds_tree->nodeChildExists(bounds_node, i)) {
          const NODE* value_child = NULL;
          if (value_node) {
            if (!value_tree->nodeHasChildren(value_node)) {
              // The value node is a leaf, just keep it the same as we
              // descend, as all we need is its value to set.
              value_child = value_node;
            } else if (value_tree->nodeChildExists(value_node, i)) {
              value_child = value_tree->getNodeChild(value_node, i);
            }
          }
          // If the value child node is NULL, there is nothing left to do,
          // unless delete_first is set.
          if (delete_first || value_child != NULL) {
            bool created_node = false;
            if (!this->nodeChildExists(node, i)) {
              // The child does not exist. Create a child if appropriate.
              if (delete_first && value_child == NULL) {
                // There is no point in creating a node that will certainly
                // be deleted. Do nothing.
              } else {
                this->createNodeChild(node, i);
                created_node = true;
              }
            }
            // It is possible to still have no child in the case that we are
            // set to delete_first and the value tree is empty here. In such a
            // case, there is no point in creating nodes just to delete them.
            if (this->nodeChildExists(node, i)) {
              key_type center_offset_key = computeCenterOffsetKey(depth, this->tree_max_val);
              OcTreeKey child_key;
              computeChildKey(i, center_offset_key, key, child_key);
              NODE* rv = setTreeValuesRecurs(this->getNodeChild(node, i),
                                             created_node,
                                             child_key,
                                             depth + 1,
                                             value_tree,
                                             bounds_tree,
                                             value_child,
                                             bounds_tree->getNodeChild(bounds_node, i),
                                             maximum_only,
                                             delete_first,
                                             copy_value_function);
              this->setNodeChild(node, i, rv);
            }
          }
        }
      }
      if (node != NULL && !this->nodeHasChildren(node)) {
        // We didn't end up with any children. Since we haven't pruned this
        // node yet, this means we lost any children we had to start with.
        // Delete.
        this->deleteNodeRecurs(node);
        node = NULL;
      }
    } else {
      // we are inside the bounds of the update
      if (delete_first) {
        if (value_node == NULL) {
          // There are no values to copy over, delete this node completely.
          this->deleteNodeRecurs(node);
          node = NULL;
        } else {
          // There is at least some value to copy over.
          // Delete all our children so they can be replaced by the values.
          this->deleteNodeChildren(node);
          // Be sure to set node_just_created as this indicates that the given
          // node is not a pruned leaf, even though it has no children.
          node_just_created = true;
        }
      }
      setTreeValuesRecurs(node, node_just_created, key, depth, value_tree, value_node, maximum_only, copy_value_function);
    }
    if (node != NULL && node_just_created && !this->nodeHasChildren(node) && value_node == NULL) {
      // We are not a leaf node and didn't end up with any children after being newly created. Delete.
      this->deleteNodeRecurs(node);
      node = NULL;
    }
    if (node != NULL) {
      // If our node is still alive, and has children, attempt a prune, and
      // then update our occupancy if we didn't prune our kids.
      if (this->nodeHasChildren(node)) {
        if (!this->pruneNode(node)) {
          node->updateOccupancyChildren();
        }
      }
    }
    return node;
  }

  template <class NODE>
  void OccupancyOcTreeBase<NODE>::setTreeValuesRecurs(NODE* node,
                                                      bool node_just_created,
                                                      const OcTreeKey& key,
                                                      unsigned int depth,
                                                      const OccupancyOcTreeBase<NODE>* value_tree,
                                                      const NODE* value_node,
                                                      bool maximum_only,
                                                      const CopyValueFunction& copy_value_function){
    if (node == NULL || value_tree == NULL || value_node == NULL)
      return;
    key_type center_offset_key = computeCenterOffsetKey(depth, this->tree_max_val);
    OcTreeKey child_key;
    if (!value_tree->nodeHasChildren(value_node)) {
      // The value node has no children (its a leaf)
      if (!this->nodeHasChildren(node)) {
        // both nodes are leafs
        if (node_just_created || !maximum_only || value_node->getValue() > node->getValue()) {
          // Update our node if it was just created (which means its empty),
          // or if we are not setting only maximum values, or if the value
          // node is a higher log odds than this node. Otherwise, there is
          // nothing to set.
          if (copy_value_function) {
            copy_value_function(value_node, node, node_just_created, key, depth);
          } else {
            node->copyData(*value_node);
          }
        }
      } else {
        // The value node is a leaf, but we are not.
        if (!maximum_only && !copy_value_function) {
          // We are updating all values, not just maximum. Delete our children
          // and copy the value node's data into our node.
          this->deleteNodeChildren(node);
          node->copyData(*value_node);
        } else {
          // Descend our tree, but keep using the node pointer from the
          // value tree's leaf to get its log odds value.  It is necessary
          // to descend. Consider, we may have log odds -1.0 and 2.0 in
          // our sub tree, and the value node may point to a leaf with
          // value of 1.0. Since we are either using maximum only, or a custom
          // copy value function, we may only change some of the child log odds
          // values.
          for (unsigned int i=0; i<8; ++i) {
            if (this->nodeChildExists(node, i)) {
              // In this case the value node is already at a leaf.
              // This value node leaf is OK to pass down the recursion, as all
              // we need is its value.
              computeChildKey(i, center_offset_key, key, child_key);
              setTreeValuesRecurs(this->getNodeChild(node, i),
                                  false,
                                  child_key,
                                  depth + 1,
                                  value_tree,
                                  value_node,
                                  true,
                                  copy_value_function);
            } else {
              // The value node has a value for this space, but we don't have
              // a child for it. Create a leaf node for our child and assign
              // it the value.
              this->createNodeChild(node, i);
              NODE* child_node = this->getNodeChild(node, i);
              if (copy_value_function) {
                // Just created the node, so set node_just_created to true
                copy_value_function(value_node, child_node, true, key, depth);
              } else {
                child_node->copyData(*value_node);
              }
            }
          }
        }
      }
    } else {
      // The value node has children. We must descend.
      for (unsigned int i=0; i<8; ++i) {
        if (value_tree->nodeChildExists(value_node, i)) {
          bool created_node = false;
          if (!this->nodeChildExists(node, i)) {
            if (!this->nodeHasChildren(node) && !node_just_created) {
              // Our node is a leaf. Expand it so we can update indivdual
              // leaves as we go.
              this->expandNode(node);
            } else {
              this->createNodeChild(node, i);
              created_node = true;
            }
          }
          computeChildKey(i, center_offset_key, key, child_key);
          setTreeValuesRecurs(this->getNodeChild(node, i),
                              created_node,
                              child_key,
                              depth + 1,
                              value_tree,
                              value_tree->getNodeChild(value_node, i),
                              maximum_only,
                              copy_value_function);
        }
      }
    }
    if (this->nodeHasChildren(node)) {
      if (!this->pruneNode(node)) {
        node->updateOccupancyChildren();
      }
    }
  }

  template <class NODE>
  void OccupancyOcTreeBase<NODE>::updateInnerOccupancy(){
    if (this->root)
      this->updateInnerOccupancyRecurs(this->root, 0);
  }

  template <class NODE>
  void OccupancyOcTreeBase<NODE>::updateInnerOccupancyRecurs(NODE* node, unsigned int depth){
    assert(node);

    // only recurse and update for inner nodes:
    if (this->nodeHasChildren(node)){
      // return early for last level:
      if (depth < this->tree_depth){
        for (unsigned int i=0; i<8; i++) {
          if (this->nodeChildExists(node, i)) {
            updateInnerOccupancyRecurs(this->getNodeChild(node, i), depth+1);
          }
        }
      }
      node->updateOccupancyChildren();
    }
  }

  template <class NODE>
  void OccupancyOcTreeBase<NODE>::toMaxLikelihood() {
    if (this->root == NULL)
      return;

    // convert bottom up
    for (unsigned int depth=this->tree_depth; depth>0; depth--) {
      toMaxLikelihoodRecurs(this->root, 0, depth);
    }

    // convert root
    nodeToMaxLikelihood(this->root);
  }

  template <class NODE>
  void OccupancyOcTreeBase<NODE>::toMaxLikelihoodRecurs(NODE* node, unsigned int depth,
      unsigned int max_depth) {

    assert(node);

    if (depth < max_depth) {
      for (unsigned int i=0; i<8; i++) {
        if (this->nodeChildExists(node, i)) {
          toMaxLikelihoodRecurs(this->getNodeChild(node, i), depth+1, max_depth);
        }
      }
    }

    else { // max level reached
      nodeToMaxLikelihood(node);
    }
  }

  template <class NODE>
  bool OccupancyOcTreeBase<NODE>::getNormals(const point3d& point, std::vector<point3d>& normals,
                                             bool unknownStatus) const {
    normals.clear();

    OcTreeKey init_key;
    if ( !OcTreeBaseImpl<NODE,AbstractOccupancyOcTree>::coordToKeyChecked(point, init_key) ) {
      OCTOMAP_WARNING_STR("Voxel out of bounds");
      return false;
    }

    // OCTOMAP_WARNING("Normal for %f, %f, %f\n", point.x(), point.y(), point.z());

    int vertex_values[8];

    OcTreeKey current_key;
    NODE* current_node;

    // There is 8 neighbouring sets
    // The current cube can be at any of the 8 vertex
    int x_index[4][4] = {{1, 1, 0, 0}, {1, 1, 0, 0}, {0, 0 -1, -1}, {0, 0 -1, -1}};
    int y_index[4][4] = {{1, 0, 0, 1}, {0, -1, -1, 0}, {0, -1, -1, 0}, {1, 0, 0, 1}};
    int z_index[2][2] = {{0, 1}, {-1, 0}};

    // Iterate over the 8 neighboring sets
    for(int m = 0; m < 2; ++m){
      for(int l = 0; l < 4; ++l){

        int k = 0;
        // Iterate over the cubes
        for(int j = 0; j < 2; ++j){
          for(int i = 0; i < 4; ++i){
            current_key[0] = init_key[0] + x_index[l][i];
            current_key[1] = init_key[1] + y_index[l][i];
            current_key[2] = init_key[2] + z_index[m][j];
            current_node = this->search(current_key);

            if(current_node){
              vertex_values[k] = this->isNodeOccupied(current_node);

              // point3d coord = this->keyToCoord(current_key);
              // OCTOMAP_WARNING_STR("vertex " << k << " at " << coord << "; value " << vertex_values[k]);
            }else{
              // Occupancy of unknown cells
              vertex_values[k] = unknownStatus;
            }
            ++k;
          }
        }

        int cube_index = 0;
        if (vertex_values[0]) cube_index |= 1;
        if (vertex_values[1]) cube_index |= 2;
        if (vertex_values[2]) cube_index |= 4;
        if (vertex_values[3]) cube_index |= 8;
        if (vertex_values[4]) cube_index |= 16;
        if (vertex_values[5]) cube_index |= 32;
        if (vertex_values[6]) cube_index |= 64;
        if (vertex_values[7]) cube_index |= 128;

        // OCTOMAP_WARNING_STR("cubde_index: " << cube_index);

        // All vertices are occupied or free resulting in no normal
        if (edgeTable[cube_index] == 0)
          return true;

        // No interpolation is done yet, we use vertexList in <MCTables.h>.
        for(int i = 0; triTable[cube_index][i] != -1; i += 3){
          point3d p1 = vertexList[triTable[cube_index][i  ]];
          point3d p2 = vertexList[triTable[cube_index][i+1]];
          point3d p3 = vertexList[triTable[cube_index][i+2]];
          point3d v1 = p2 - p1;
          point3d v2 = p3 - p1;

          // OCTOMAP_WARNING("Vertex p1 %f, %f, %f\n", p1.x(), p1.y(), p1.z());
          // OCTOMAP_WARNING("Vertex p2 %f, %f, %f\n", p2.x(), p2.y(), p2.z());
          // OCTOMAP_WARNING("Vertex p3 %f, %f, %f\n", p3.x(), p3.y(), p3.z());

          // Right hand side cross product to retrieve the normal in the good
          // direction (pointing to the free nodes).
          normals.push_back(v1.cross(v2).normalize());
        }
      }
    }

    return true;
  }

  template <class NODE>
  bool OccupancyOcTreeBase<NODE>::castRay(const point3d& origin, const point3d& directionP, point3d& end,
                                          bool ignoreUnknown, double maxRange) const {

    /// ----------  see OcTreeBase::computeRayKeys  -----------

    // Initialization phase -------------------------------------------------------
    OcTreeKey current_key;
    if ( !OcTreeBaseImpl<NODE,AbstractOccupancyOcTree>::coordToKeyChecked(origin, current_key) ) {
      OCTOMAP_WARNING_STR("Coordinates out of bounds during ray casting");
      return false;
    }

    NODE* startingNode = this->search(current_key);
    if (startingNode){
      if (this->isNodeOccupied(startingNode)){
        // Occupied node found at origin
        // (need to convert from key, since origin does not need to be a voxel center)
        end = this->keyToCoord(current_key);
        return true;
      }
    } else if(!ignoreUnknown){
      end = this->keyToCoord(current_key);
      return false;
    }

    point3d direction = directionP.normalized();
    bool max_range_set = (maxRange > 0.0);

    int step[3];
    double tMax[3];
    double tDelta[3];

    for(unsigned int i=0; i < 3; ++i) {
      // compute step direction
      if (direction(i) > 0.0) step[i] =  1;
      else if (direction(i) < 0.0)   step[i] = -1;
      else step[i] = 0;

      // compute tMax, tDelta
      if (step[i] != 0) {
        // corner point of voxel (in direction of ray)
        double voxelBorder = this->keyToCoord(current_key[i]);
        voxelBorder += double(step[i] * this->resolution * 0.5);

        tMax[i] = ( voxelBorder - origin(i) ) / direction(i);
        tDelta[i] = this->resolution / fabs( direction(i) );
      }
      else {
        tMax[i] =  std::numeric_limits<double>::max();
        tDelta[i] = std::numeric_limits<double>::max();
      }
    }

    if (step[0] == 0 && step[1] == 0 && step[2] == 0){
    	OCTOMAP_ERROR("Raycasting in direction (0,0,0) is not possible!");
    	return false;
    }

    // for speedup:
    double maxrange_sq = maxRange *maxRange;

    // Incremental phase  ---------------------------------------------------------

    bool done = false;

    while (!done) {
      unsigned int dim;

      // find minimum tMax:
      if (tMax[0] < tMax[1]){
        if (tMax[0] < tMax[2]) dim = 0;
        else                   dim = 2;
      }
      else {
        if (tMax[1] < tMax[2]) dim = 1;
        else                   dim = 2;
      }

      // check for overflow:
      if ((step[dim] < 0 && current_key[dim] == 0)
    		  || (step[dim] > 0 && current_key[dim] == 2* this->tree_max_val-1))
      {
        OCTOMAP_WARNING("Coordinate hit bounds in dim %d, aborting raycast\n", dim);
        // return border point nevertheless:
        end = this->keyToCoord(current_key);
        return false;
      }

      // advance in direction "dim"
      current_key[dim] += step[dim];
      tMax[dim] += tDelta[dim];


      // generate world coords from key
      end = this->keyToCoord(current_key);

      // check for maxrange:
      if (max_range_set){
        double dist_from_origin_sq(0.0);
        for (unsigned int j = 0; j < 3; j++) {
          dist_from_origin_sq += ((end(j) - origin(j)) * (end(j) - origin(j)));
        }
        if (dist_from_origin_sq > maxrange_sq)
          return false;

      }

      NODE* currentNode = this->search(current_key);
      if (currentNode){
        if (this->isNodeOccupied(currentNode)) {
          done = true;
          break;
        }
        // otherwise: node is free and valid, raycasting continues
      } else if (!ignoreUnknown){ // no node found, this usually means we are in "unknown" areas
        return false;
      }
    } // end while

    return true;
  }

  template <class NODE>
  bool OccupancyOcTreeBase<NODE>::getRayIntersection (const point3d& origin, const point3d& direction, const point3d& center,
                 point3d& intersection, double delta/*=0.0*/) const {
    // We only need three normals for the six planes
    octomap::point3d normalX(1, 0, 0);
    octomap::point3d normalY(0, 1, 0);
    octomap::point3d normalZ(0, 0, 1);

    // One point on each plane, let them be the center for simplicity
    octomap::point3d pointXNeg(center(0) - double(this->resolution / 2.0), center(1), center(2));
    octomap::point3d pointXPos(center(0) + double(this->resolution / 2.0), center(1), center(2));
    octomap::point3d pointYNeg(center(0), center(1) - double(this->resolution / 2.0), center(2));
    octomap::point3d pointYPos(center(0), center(1) + double(this->resolution / 2.0), center(2));
    octomap::point3d pointZNeg(center(0), center(1), center(2) - double(this->resolution / 2.0));
    octomap::point3d pointZPos(center(0), center(1), center(2) + double(this->resolution / 2.0));

    double lineDotNormal = 0.0;
    double d = 0.0;
    double outD = std::numeric_limits<double>::max();
    octomap::point3d intersect;
    bool found = false;

    // Find the intersection (if any) with each place
    // Line dot normal will be zero if they are parallel, in which case no intersection can be the entry one
    // if there is an intersection does it occur in the bounded plane of the voxel
    // if yes keep only the closest (smallest distance to sensor origin).
    if((lineDotNormal = normalX.dot(direction)) != 0.0){   // Ensure lineDotNormal is non-zero (assign and test)
      d = (pointXNeg - origin).dot(normalX) / lineDotNormal;
      intersect = direction * double(d) + origin;
      if(!(intersect(1) < (pointYNeg(1) - 1e-6) || intersect(1) > (pointYPos(1) + 1e-6) ||
         intersect(2) < (pointZNeg(2) - 1e-6) || intersect(2) > (pointZPos(2) + 1e-6))){
        outD = std::min(outD, d);
        found = true;
      }

      d = (pointXPos - origin).dot(normalX) / lineDotNormal;
      intersect = direction * double(d) + origin;
      if(!(intersect(1) < (pointYNeg(1) - 1e-6) || intersect(1) > (pointYPos(1) + 1e-6) ||
         intersect(2) < (pointZNeg(2) - 1e-6) || intersect(2) > (pointZPos(2) + 1e-6))){
        outD = std::min(outD, d);
        found = true;
      }
    }

    if((lineDotNormal = normalY.dot(direction)) != 0.0){   // Ensure lineDotNormal is non-zero (assign and test)
      d = (pointYNeg - origin).dot(normalY) / lineDotNormal;
      intersect = direction * double(d) + origin;
      if(!(intersect(0) < (pointXNeg(0) - 1e-6) || intersect(0) > (pointXPos(0) + 1e-6) ||
         intersect(2) < (pointZNeg(2) - 1e-6) || intersect(2) > (pointZPos(2) + 1e-6))){
        outD = std::min(outD, d);
        found = true;
      }

      d = (pointYPos - origin).dot(normalY) / lineDotNormal;
      intersect = direction * double(d) + origin;
      if(!(intersect(0) < (pointXNeg(0) - 1e-6) || intersect(0) > (pointXPos(0) + 1e-6) ||
         intersect(2) < (pointZNeg(2) - 1e-6) || intersect(2) > (pointZPos(2) + 1e-6))){
        outD = std::min(outD, d);
        found = true;
      }
    }

    if((lineDotNormal = normalZ.dot(direction)) != 0.0){   // Ensure lineDotNormal is non-zero (assign and test)
      d = (pointZNeg - origin).dot(normalZ) / lineDotNormal;
      intersect = direction * double(d) + origin;
      if(!(intersect(0) < (pointXNeg(0) - 1e-6) || intersect(0) > (pointXPos(0) + 1e-6) ||
         intersect(1) < (pointYNeg(1) - 1e-6) || intersect(1) > (pointYPos(1) + 1e-6))){
        outD = std::min(outD, d);
        found = true;
      }

      d = (pointZPos - origin).dot(normalZ) / lineDotNormal;
      intersect = direction * double(d) + origin;
      if(!(intersect(0) < (pointXNeg(0) - 1e-6) || intersect(0) > (pointXPos(0) + 1e-6) ||
         intersect(1) < (pointYNeg(1) - 1e-6) || intersect(1) > (pointYPos(1) + 1e-6))){
        outD = std::min(outD, d);
        found = true;
      }
    }

    // Substract (add) a fraction to ensure no ambiguity on the starting voxel
    // Don't start on a boundary.
    if(found)
      intersection = direction * double(outD + delta) + origin;

    return found;
  }


  template <class NODE> inline bool
  OccupancyOcTreeBase<NODE>::integrateMissOnRay(const point3d& origin, const point3d& end, bool lazy_eval) {

    if (!this->computeRayKeys(origin, end, this->keyrays.at(0))) {
      return false;
    }

    for(KeyRay::iterator it=this->keyrays[0].begin(); it != this->keyrays[0].end(); it++) {
      updateNode(*it, false, lazy_eval); // insert freespace measurement
    }

    return true;
  }

  template <class NODE> bool
  OccupancyOcTreeBase<NODE>::insertRay(const point3d& origin, const point3d& end, double maxrange, bool lazy_eval)
  {
    // cut ray at maxrange
    if ((maxrange > 0) && ((end - origin).norm () > maxrange))
      {
        point3d direction = (end - origin).normalized ();
        point3d new_end = origin + direction * maxrange;
        return integrateMissOnRay(origin, new_end,lazy_eval);
      }
    // insert complete ray
    else
      {
        if (!integrateMissOnRay(origin, end,lazy_eval))
          return false;
        updateNode(end, true, lazy_eval); // insert hit cell
        return true;
      }
  }

  template <class NODE>
  void OccupancyOcTreeBase<NODE>::setBBXMin (point3d& min) {
    bbx_min = min;
    if (!this->coordToKeyChecked(bbx_min, bbx_min_key)) {
      OCTOMAP_ERROR("ERROR while generating bbx min key.\n");
    }
  }

  template <class NODE>
  void OccupancyOcTreeBase<NODE>::setBBXMax (point3d& max) {
    bbx_max = max;
    if (!this->coordToKeyChecked(bbx_max, bbx_max_key)) {
      OCTOMAP_ERROR("ERROR while generating bbx max key.\n");
    }
  }


  template <class NODE>
  bool OccupancyOcTreeBase<NODE>::inBBX(const point3d& p) const {
    return ((p.x() >= bbx_min.x()) && (p.y() >= bbx_min.y()) && (p.z() >= bbx_min.z()) &&
            (p.x() <= bbx_max.x()) && (p.y() <= bbx_max.y()) && (p.z() <= bbx_max.z()) );
  }


  template <class NODE>
  bool OccupancyOcTreeBase<NODE>::inBBX(const OcTreeKey& key) const {
    return ((key[0] >= bbx_min_key[0]) && (key[1] >= bbx_min_key[1]) && (key[2] >= bbx_min_key[2]) &&
            (key[0] <= bbx_max_key[0]) && (key[1] <= bbx_max_key[1]) && (key[2] <= bbx_max_key[2]) );
  }

  template <class NODE>
  point3d OccupancyOcTreeBase<NODE>::getBBXBounds () const {
    octomap::point3d obj_bounds = (bbx_max - bbx_min);
    obj_bounds /= 2.;
    return obj_bounds;
  }

  template <class NODE>
  point3d OccupancyOcTreeBase<NODE>::getBBXCenter () const {
    octomap::point3d obj_bounds = (bbx_max - bbx_min);
    obj_bounds /= 2.;
    return bbx_min + obj_bounds;
  }

  // -- I/O  -----------------------------------------

  template <class NODE>
  std::istream& OccupancyOcTreeBase<NODE>::readBinaryData(std::istream &s){
    // tree needs to be newly created or cleared externally
    if (this->root) {
      OCTOMAP_ERROR_STR("Trying to read into an existing tree.");
      return s;
    }

    this->root = this->allocNode();
    this->readBinaryNode(s, this->root);
    this->size_changed = true;
    this->tree_size = OcTreeBaseImpl<NODE,AbstractOccupancyOcTree>::calcNumNodes();  // compute number of nodes
    return s;
  }

  template <class NODE>
  std::ostream& OccupancyOcTreeBase<NODE>::writeBinaryData(std::ostream &s) const{
    OCTOMAP_DEBUG("Writing %zu nodes to output stream...", this->size());
    if (this->root)
      this->writeBinaryNode(s, this->root);
    return s;
  }

  template <class NODE>
  std::istream& OccupancyOcTreeBase<NODE>::readBinaryNode(std::istream &s, NODE* node){

    assert(node);

    char child1to4_char;
    char child5to8_char;
    s.read((char*)&child1to4_char, sizeof(char));
    s.read((char*)&child5to8_char, sizeof(char));

    std::bitset<8> child1to4 ((unsigned long long) child1to4_char);
    std::bitset<8> child5to8 ((unsigned long long) child5to8_char);

    //     std::cout << "read:  "
    //        << child1to4.to_string<char,std::char_traits<char>,std::allocator<char> >() << " "
    //        << child5to8.to_string<char,std::char_traits<char>,std::allocator<char> >() << std::endl;


    // inner nodes default to occupied
    node->setLogOdds(this->clamping_thres_max);

    for (unsigned int i=0; i<4; i++) {
      if ((child1to4[i*2] == 1) && (child1to4[i*2+1] == 0)) {
        // child is free leaf
        this->createNodeChild(node, i);
        this->getNodeChild(node, i)->setLogOdds(this->clamping_thres_min);
      }
      else if ((child1to4[i*2] == 0) && (child1to4[i*2+1] == 1)) {
        // child is occupied leaf
        this->createNodeChild(node, i);
        this->getNodeChild(node, i)->setLogOdds(this->clamping_thres_max);
      }
      else if ((child1to4[i*2] == 1) && (child1to4[i*2+1] == 1)) {
        // child has children
        this->createNodeChild(node, i);
        this->getNodeChild(node, i)->setLogOdds(-200.); // child is unkown, we leave it uninitialized
      }
    }
    for (unsigned int i=0; i<4; i++) {
      if ((child5to8[i*2] == 1) && (child5to8[i*2+1] == 0)) {
        // child is free leaf
        this->createNodeChild(node, i+4);
        this->getNodeChild(node, i+4)->setLogOdds(this->clamping_thres_min);
      }
      else if ((child5to8[i*2] == 0) && (child5to8[i*2+1] == 1)) {
        // child is occupied leaf
        this->createNodeChild(node, i+4);
        this->getNodeChild(node, i+4)->setLogOdds(this->clamping_thres_max);
      }
      else if ((child5to8[i*2] == 1) && (child5to8[i*2+1] == 1)) {
        // child has children
        this->createNodeChild(node, i+4);
        this->getNodeChild(node, i+4)->setLogOdds(-200.); // set occupancy when all children have been read
      }
      // child is unkown, we leave it uninitialized
    }

    // read children's children and set the label
    for (unsigned int i=0; i<8; i++) {
      if (this->nodeChildExists(node, i)) {
        NODE* child = this->getNodeChild(node, i);
        if (fabs(child->getLogOdds() + 200.)<1e-3) {
          readBinaryNode(s, child);
          child->setLogOdds(child->getMaxChildLogOdds());
        }
      } // end if child exists
    } // end for children

    return s;
  }

  template <class NODE>
  std::ostream& OccupancyOcTreeBase<NODE>::writeBinaryNode(std::ostream &s, const NODE* node) const{

    assert(node);

    // 2 bits for each children, 8 children per node -> 16 bits
    std::bitset<8> child1to4;
    std::bitset<8> child5to8;

    // 10 : child is free node
    // 01 : child is occupied node
    // 00 : child is unkown node
    // 11 : child has children


    // speedup: only set bits to 1, rest is init with 0 anyway,
    //          can be one logic expression per bit

    for (unsigned int i=0; i<4; i++) {
      if (this->nodeChildExists(node, i)) {
        const NODE* child = this->getNodeChild(node, i);
        if      (this->nodeHasChildren(child))  { child1to4[i*2] = 1; child1to4[i*2+1] = 1; }
        else if (this->isNodeOccupied(child)) { child1to4[i*2] = 0; child1to4[i*2+1] = 1; }
        else                            { child1to4[i*2] = 1; child1to4[i*2+1] = 0; }
      }
      else {
        child1to4[i*2] = 0; child1to4[i*2+1] = 0;
      }
    }

    for (unsigned int i=0; i<4; i++) {
      if (this->nodeChildExists(node, i+4)) {
        const NODE* child = this->getNodeChild(node, i+4);
        if      (this->nodeHasChildren(child))  { child5to8[i*2] = 1; child5to8[i*2+1] = 1; }
        else if (this->isNodeOccupied(child)) { child5to8[i*2] = 0; child5to8[i*2+1] = 1; }
        else                            { child5to8[i*2] = 1; child5to8[i*2+1] = 0; }
      }
      else {
        child5to8[i*2] = 0; child5to8[i*2+1] = 0;
      }
    }
    //     std::cout << "wrote: "
    //        << child1to4.to_string<char,std::char_traits<char>,std::allocator<char> >() << " "
    //        << child5to8.to_string<char,std::char_traits<char>,std::allocator<char> >() << std::endl;

    char child1to4_char = (char) child1to4.to_ulong();
    char child5to8_char = (char) child5to8.to_ulong();

    s.write((char*)&child1to4_char, sizeof(char));
    s.write((char*)&child5to8_char, sizeof(char));

    // write children's children
    for (unsigned int i=0; i<8; i++) {
      if (this->nodeChildExists(node, i)) {
        const NODE* child = this->getNodeChild(node, i);
        if (this->nodeHasChildren(child)) {
          writeBinaryNode(s, child);
        }
      }
    }

    return s;
  }

  //-- Occupancy queries on nodes:

  template <class NODE>
  void OccupancyOcTreeBase<NODE>::updateNodeLogOdds(NODE* occupancyNode, const float& update) const {
    occupancyNode->addValue(update);
    if (occupancyNode->getLogOdds() < this->clamping_thres_min) {
      occupancyNode->setLogOdds(this->clamping_thres_min);
      return;
    }
    if (occupancyNode->getLogOdds() > this->clamping_thres_max) {
      occupancyNode->setLogOdds(this->clamping_thres_max);
    }
  }

  template <class NODE>
  void OccupancyOcTreeBase<NODE>::updateNodeLogOddsAndTrackChanges(NODE* node, const float& update,
                                                                   bool node_just_created, const OcTreeKey& key,
                                                                   unsigned int depth) {
    if (use_change_detection || nodeValueChangeCallback) {
      bool occBefore = this->isNodeOccupied(node);
      float valBefore = node->getLogOdds();
      updateNodeLogOdds(node, update);

      if(node->getLogOdds() != valBefore || node_just_created)
      {
        valueChangeCallbackWrapper(key, depth, node_just_created, valBefore, occBefore, node->getLogOdds(), this->isNodeOccupied(node));
      }

      if (use_change_detection) {
        if (node_just_created){  // new node
          changed_keys.insert(std::pair<OcTreeKey,bool>(key, true));
        } else if (occBefore != this->isNodeOccupied(node)) {  // occupancy changed, track it
          KeyBoolMap::iterator it = changed_keys.find(key);
          if (it == changed_keys.end())
            changed_keys.insert(std::pair<OcTreeKey,bool>(key, false));
          else if (it->second == false)
            changed_keys.erase(it);
        }
      }
    } else {
      updateNodeLogOdds(node, update);
    }
  }

  template <class NODE>
  void OccupancyOcTreeBase<NODE>::integrateHit(NODE* occupancyNode) const {
    updateNodeLogOdds(occupancyNode, this->prob_hit_log);
  }

  template <class NODE>
  void OccupancyOcTreeBase<NODE>::integrateMiss(NODE* occupancyNode) const {
    updateNodeLogOdds(occupancyNode, this->prob_miss_log);
  }

  template <class NODE>
  void OccupancyOcTreeBase<NODE>::nodeToMaxLikelihood(NODE* occupancyNode) const{
    if (this->isNodeOccupied(occupancyNode))
      occupancyNode->setLogOdds(this->clamping_thres_max);
    else
      occupancyNode->setLogOdds(this->clamping_thres_min);
  }

  template <class NODE>
  void OccupancyOcTreeBase<NODE>::nodeToMaxLikelihood(NODE& occupancyNode) const{
    if (this->isNodeOccupied(occupancyNode))
      occupancyNode.setLogOdds(this->clamping_thres_max);
    else
      occupancyNode.setLogOdds(this->clamping_thres_min);
  }

  template <class NODE>
  void OccupancyOcTreeBase<NODE>::valueChangeCallbackWrapper(const OcTreeKey& key, unsigned int depth, const bool node_just_created,
      const float prev_full_val, const bool prev_binary_val,
      const float curr_full_val, const bool curr_binary_val) {
    if (nodeValueChangeCallback)
      nodeValueChangeCallback(key, depth, node_just_created, prev_full_val, prev_binary_val, curr_full_val, curr_binary_val);
  }

} // namespace
